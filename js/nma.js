import { buildKirchhoff } from './geometry.js';
import { F32toHeap } from './wasm-utils.js';


/**
 * Diagonalizes a matrix using WASM-compiled code.
 * @param {Float32Array} mtx - (full) matrix in flat form.
 * @param {number} d - dimension of the (square) matrix.
 * @param {number} m - number of modes to calculate/return.
 */
export function diagonalize(mtx, d, m) {

    let _m = m + 1;  // account for 0th mode

    let _evals = new Float32Array(_m);
    let _evecs = new Float32Array(_m * d);

    let _kmtx = F32toHeap(mtx);
    let __evals = F32toHeap(_evals);
    let __evecs = F32toHeap(_evecs);

    console.time('diag')
    Module._diagonalize(_kmtx, d, _m, __evals, __evecs)

    // Get data from WASM heap
    // Eigenvalues (ignore 0th mode - eigval = 0)
    let evals = new Float32Array(m);
    for (let i = 1; i < _m; i++) {
        let ieval = Module.HEAPF32[__evals / Float32Array.BYTES_PER_ELEMENT + i]

        if (ieval < 1e-5) {
            throw("Unexpected zero eigenvalue for mode " + i + ": " + ieval);
        }

        evals[i - 1] = ieval
    }

    // Get full eigenvector (d entries per each m mode)
    // Array is in column-major order
    let evecs = []
    for (let i = d; i < _m*d; i++) {  // start at d to ignore first mode
        let idx = (i % d) * _m + Math.floor(i / d);
        evecs.push(Module.HEAPF32[__evecs / Float32Array.BYTES_PER_ELEMENT + idx])
    }

    // Now free memory
    Module._free(_kmtx);
    Module._free(_evals);
    Module._free(_evecs);

    console.log('eigenvalues:', evals);
    console.log('eigenvectors:', evecs);
    console.timeEnd('diag')

    return [evals, evecs]
}

export function GNM(model, nmodes) {

    // Get flattened coordinate array from model
    let coordinates = [];
    let selection = [];
    let resid2bead = {};
    for (let b = 0; b < model.length; b++) {

        let resid = model[b].resi + model[b].chain;

        if (model[b].atom == "CA" || model[b].atom == "C4*") {
            coordinates.push(model[b].x)
            coordinates.push(model[b].y)
            coordinates.push(model[b].z)

            selection.push(resid);
        }

        if (!(resid in resid2bead)) {
            resid2bead[resid] = []
        }
        resid2bead[resid].push(b)
    }
    coordinates = Float32Array.from(coordinates);

    const [natoms, kmtx] = buildKirchhoff(coordinates);
    const [evals, evecs] = diagonalize(kmtx, natoms, nmodes);

    // Get square fluctuations
    const sqfluctuations = new Float32Array(natoms);
    const variance = evals.map(ev => (1/ev));
    for (let i = 0; i < evals.length; i++) {
        for (let a = 0; a < natoms; a++) {
            sqfluctuations[a] += (evecs[i * natoms + a] * evecs[i * natoms + a]) * variance[i]
        }
    }

    // Assign fluctuations to b factor of model
    for (let a = 0; a < natoms; a++) {
        // Propagate to atoms of same residue.
        let sele = selection[a];
        let atoms = resid2bead[sele];
        for (let i = 0; i < atoms.length; i++) {
            let idx = atoms[i];
            model[idx].b = sqfluctuations[a]
        }
    }

    return sqfluctuations;
}
