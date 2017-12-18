// Functions to perform GNM on a structure

// Constants
const sq_threshold = 100; // squared distance threshold for Kirchoff Matrix
const dtype = jsfeat.F32_t | jsfeat.C1_t; // data type for jsfeat matrices
const selectedModes = [0, 1, 2, 3, 4, 5];

// Global variables
let mtxOriginal; // Original Kirchoff Matrix to keep track of connectivity
let mtx; // Kirchoff Matrix
let numAtoms; // Number of atoms in the structure

let eigenvalues;
let eigenvectors;
let sqfluctuations;
let overallFluctuation;
let initialMaximumFluctuation;


let residueToIdx = {}; // residue number to mtx index

// Wrapper function called from main.js
function NMA(mtx) {
    GNM(mtx);
    calcSquareFluctuations(selectedModes);
}

// Auxiliary Functions

// Iterates over the CA/C5' of a structure
// and calculates the connectivity (Kirchhoff) matrix
function buildKirchhoff(structure) {

    const atoms = [];
    const atomsel = new NGL.Selection(".CA or .C5'")
    component.structure.eachAtom(function(a) {
        atoms.push([a.x, a.y, a.z]);
        residueToIdx[a.resno] = atoms.length - 1;
    }, atomsel);

    numAtoms = atoms.length;
    console.log("Selected " + numAtoms + " atoms");

    // Allocate Matrix
    mtx = new jsfeat.matrix_t(numAtoms, numAtoms, dtype);

    // Do distance calculation
    let t0 = performance.now();
    
    let d_ij_sq;
    for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
            d_ij_sq = sqEuclideanDistance(atoms[i], atoms[j]);
            // console.log(i, j, d_ij_sq);
            if (d_ij_sq <= sq_threshold) {
                mtx.data[i*numAtoms + j] = -1;
                mtx.data[i*numAtoms + i] += 1;
                mtx.data[j*numAtoms + i] = -1;
                mtx.data[j*numAtoms + j] += 1;
            }
        }
    }

    t1 = performance.now();
    console.log("Calculated distances in: " + (t1 - t0) + " milliseconds");

    // Store initial matrix (not very efficient, but this is the 21st century. RAM is cheap)
    mtxOriginal = new jsfeat.matrix_t(numAtoms, numAtoms, dtype);
    for (let i = 0; i < mtx.data.length; i++) {
        mtxOriginal.data[i] = mtx.data[i];
    }
}


// Alters the matrix 
// e.g. idxA = 0, idxB = 1 <changes> (0,1), (1,0) (0,0) (1,1)
// If factor is negative, will restore the original value of the matrix.
function editKirchhoff(idxA, idxB, factor = 10) {

        if (factor < 0) {
            mtx.data[idxB * numAtoms + idxA] = mtxOriginal.data[idxB * numAtoms + idxA];
            mtx.data[idxA * numAtoms + idxB] = mtxOriginal.data[idxA * numAtoms + idxB];
        }
        else {
            mtx.data[idxB * numAtoms + idxA] = -1 * factor;
            mtx.data[idxA * numAtoms + idxB] = -1 * factor;            
        }

        // Sum all columns/rows to get the diagonal right
        // Slower than setting directly but easier to 'undo'
        let diagonal_A = (idxA * numAtoms + idxA);
        mtx.data[diagonal_A] = getColSum(idxA);

        let diagonal_B = (idxB * numAtoms + idxB);
        mtx.data[diagonal_B] = getColSum(idxB);
}

function getColSum(idx) {
    let sum_col = 0;
    let diagonal = (idx * numAtoms + idx);
    for (let c = 0; c < numAtoms; c++) {
        let cell = (c * numAtoms) + idx;
        if (cell === diagonal) { 
            continue
        };
        sum_col += Math.abs(mtx.data[cell]);
    }
    return sum_col;
}

// Eigendecomposition from JSFeat
// Fills the two eigenvectors/eigenvalues arrays
function GNM(kirchoffMatrix, modes = 20) {

    // let sumdiag = 0.0;
    // for (let x = 0; x < kirchoffMatrix.rows; x++) {
    //     sumdiag += kirchoffMatrix.data[(x * numAtoms + x)];
    // }
    // console.log("sum of diagonal elements: " + sumdiag);

    // Pre-allocate result arrays
    let evals = new jsfeat.matrix_t(numAtoms, 1, dtype);
    let evecs = new jsfeat.matrix_t(numAtoms, numAtoms, dtype);
    
    t0 = performance.now();
    jsfeat.linalg.eigenVV(mtx, evecs, evals); // perf. bottleneck
    t1 = performance.now();
    console.log("Diagonalized matrix in: " + (t1 - t0) + " milliseconds");
    
    // Check for zero-value modes (should be only 1)
    const small = 1e-100; // zero enough
    let numZeroModes = sumArray(evals.data.map(eval => (eval <= small) ? 1: 0));
    console.log(numZeroModes + " zero modes detected");
    if (numZeroModes === 0) {
        numZeroModes = 1;
    }

    // Modes are sorted in decreasing order of eigenvalue
    // evecs.data contains a flattened matrix:
    // [a, b, c, d, ...]
    //  |________|
    //   numAtoms = per vector
    //
    t0 = performance.now();
    eigenvectors = [];
    let mode = [];
    for (let i = 0; i <= evecs.data.length; i++) {
        if (mode.length > 0 && (i % numAtoms) == 0) {
            eigenvectors.push(mode)
            mode = [];
        }
        mode.push(evecs.data[i])
    }
    eigenvectors.reverse();
    eigenvectors = eigenvectors.slice(numZeroModes, numZeroModes + modes);

    // Sort, slice, and pass eigenvalues to Array
    eigenvalues = [].slice.call(evals.data);
    eigenvalues.reverse();
    eigenvalues = eigenvalues.slice(numZeroModes, numZeroModes + modes);
    
    t1 = performance.now();
    console.log("Sorted eigenvals/vecs in: " + (t1 - t0) + " milliseconds");
}

// Calculates atomic square fluctuations from calculated modes
// Parameter is the index/indices of the modes we want to use in the calculation
// Expects two global variables: eigenvalues/eigenvectors to be defined.
function calcSquareFluctuations(modeIndices) {

    if (eigenvectors === undefined || eigenvalues === undefined) {
        throw "Eigenvectors or Eigenvalues are not defined";
    }

    const modeVariance = eigenvalues.map(ev => (1/ev));
    sqfluctuations = new Array(numAtoms).fill(0);

    t0 = performance.now();
    if (modeIndices.length == 1) { // only one
        const m = modeIndices[0];
        sqfluctuations = entrywiseMultiplication(elementwisePow(eigenvectors[m], 2), modeVariance[m]);
    }
    else {
        for (let i = 0; i < modeIndices.length; i++) {
            let m = modeIndices[i];
            const modeFluct = entrywiseMultiplication(elementwisePow(eigenvectors[m], 2), modeVariance[m]);
            sqfluctuations = sqfluctuations.map((num, idx) => num + modeFluct[idx]);
        }
    }

    // Extend fluctuations to ALL atoms (HETATM as 0.0)
    all_sqfluctuations = [];
    let res_idx = 0;
    component.structure.eachResidue(function(r) {
        
        let value = 0.0;
        
        if (r.hetero == false) {
            value = sqfluctuations[res_idx];
            res_idx = res_idx + 1;
        }

        for (let i = 0; i < r.atomCount; i++) {
            all_sqfluctuations.push(value)  
        }
    });

    sqfluctuations = all_sqfluctuations;

    // overallFluctuation = sumArray(sqfluctuations) / sqfluctuations.length;

    // Get initial max fluct for coloring purposes
    if (initialMaximumFluctuation === undefined) {
        initialMaximumFluctuation = Math.max.apply(null, sqfluctuations);
    }


    t1 = performance.now();
    console.log("Calculated atomic square fluctuations in: " + (t1 - t0) + " milliseconds");
}

// Obvious what this does
function sqEuclideanDistance(a, b) {
    return ((a[0] - b[0])*(a[0] - b[0]) + 
            (a[1] - b[1])*(a[1] - b[1]) + 
            (a[2] - b[2])*(a[2] - b[2]))
}

// Sums all elements of an array
function sumArray(array) {
    const summer = (accumulator, i) => accumulator + i;
    return array.reduce(summer)
}

// [1, 2, 3] * 2 = [2, 4, 6]
function entrywiseMultiplication(array, scalar) {
    result = [];
    for (let i = 0; i < array.length; i++) {
        result.push(array[i] * scalar);
    }
    return result;
}

// [1, 2, 3] ^ 2 = [1, 4, 9] 
function elementwisePow(array, exponent) {
    result = [];
    for (let i = 0; i < array.length; i++) {
        result.push(Math.pow(array[i], exponent));
    }
    return result;  
}


// Debug function
function prettyMatrix() {
    let temp = [];
    for (let i = 0; i <= mtx.data.length; i++) {
        if (i > 0 && (i % mtx.rows == 0)) {
            console.log(temp)
            temp = [];
        }
        temp.push(mtx.data[i]);
    }
}
