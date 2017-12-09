//
// Define global variables
//
let stage;
let component;
let numAtoms;
let elapsedTime;

// NMA variables
let eigenvalues;
let eigenvectors;
let modeVariance;

// Regexes
const fullStructure = new RegExp("[0-9A-Za-z]{4}$");
const partialStructure = new RegExp("([0-9A-Za-z]{4})\.([A-Z0-9a-z]{1})$");

// Constants
let numModes = 10; // Can be changed later on user input
const sq_threshold = 100; // squared distance threshold for Kirchoff Matrix
//

// On DOM loaded: setup PDB and plots
document.addEventListener("DOMContentLoaded", function () {
	setupViewer();
});

// Main function executed on pressing enter in the input box (char code 13)
async function RunApp(e) {

	if (e.keyCode != 13) {
		return true;
	}
	else {
		// Clear and reveal progress div
		document.getElementById("progress_div").innerHTML = "";
		document.getElementById("progress_div").className = "visible"
		// Same for stats div
		document.getElementById("stats").innerHTML = "";
	}

	// Clear existing data
	stage.removeAllComponents();

	var structureCode = document.getElementById('input_pdb').value;
	
	// Validate PDB code and build URL
	// https://webchem.ncbr.muni.cz/CoordinateServer/1ggr/chains?modelId=1&atomSitesOnly=1&authAsymId=A
	// https://webchem.ncbr.muni.cz/CoordinateServer/1ggr/full?modelId=1&atomSitesOnly=1

	let url = "";
	if (fullStructure.test(structureCode)) {
		url = "https://webchem.ncbr.muni.cz/CoordinateServer/" + structureCode + "/full?modelId=1&atomSitesOnly=1";
	} 
	else if (partialStructure.test(structureCode)) {
		let match = partialStructure.exec(structureCode);
		url = "https://webchem.ncbr.muni.cz/CoordinateServer/" + match[1] + "/chains?modelId=1&atomSitesOnly=1&authAsymId=" + match[2];	
	}
	else {
		await updateProgress("Hmm.. that does not look like a valid PDB code: '" + structureCode + "'");
		return false
	}

	let startTime = performance.now();
	
	// Load Structure
	await updateProgress("Fetching your favorite protein");
	component = await stage.loadFile(url, { ext: "cif" });

	await updateProgress("Making it look pretty");
	await renderMolecule()  // Render molecule as hyperballs. Return structure object

	await updateProgress("Doing all the heavy work");
	await doNMA(numModes)
	elapsedTime = ((performance.now() - startTime) / 1000).toFixed(2);

	// Randomly perturb coordinates
	// setInterval(updateCoords, 1000);

	// Finally, fade out load box
	await updateProgress("Ta-dah! Have fun!", 1000);

	document.getElementById("loadbox").className = "hidden";
	document.getElementById("progress_div").className = "invisible"

	// Update stats div
	let structureNumbers = "PDB ID: " + structureCode + " | No. Atoms: " + component.structure.atomCount + " | Elapsed Time: " + elapsedTime + " seconds";
	document.getElementById("stats").innerHTML = structureNumbers;

	return true;
}

function updateProgress(message, delay) {

	let el = document.getElementById('progress_div');
	el.innerHTML = message;

	// Hackish way to allow DOM/UI updates
	return new Promise((resolve, reject) => {
		setTimeout(() => resolve(true), delay);
	});
}

function renderMolecule() {

	component.addRepresentation('hyperball', { sele: 'protein', 
											   color: 'atomindex',
											   quality: 'high' }); 

	var pa = component.structure.getView(new NGL.Selection('.CA')).getPrincipalAxes();
	stage.animationControls.rotate(pa.getRotationQuaternion(), 0);
	stage.autoView();
}

// Setup NGL stage and paramters
function setupViewer() {
	
	var stageParams = { backgroundColor: "white" };
	stage = new NGL.Stage("ngl-canvas", stageParams);

	// Handle resizing events
	function handleResize() {
		stage.handleResize();
	}

	window.addEventListener("resize", handleResize, false);
};

// NMA Calculation Code
function doNMA(numModes) {

	const atoms = []; // c-alpha array
	const atomsel = new NGL.Selection(".CA")
	component.structure.eachAtom(function(a) {
		atoms.push([a.x, a.y, a.z]);
	}, atomsel);
	// console.log("Selected " + atoms.length + " atoms");

	numAtoms = atoms.length;

	// Build Kirchoff Matrix
	let mtx = zeros(numAtoms, numAtoms);

	let t0 = performance.now();
	let d_ij_sq;
	for (let i = 0; i < atoms.length; i++) {
		for (let j = i + 1; j < atoms.length; j++) {
			d_ij_sq = sqEuclideanDistance(atoms[i], atoms[j]);
			// console.log(i, j, d_ij_sq);
			if (d_ij_sq <= sq_threshold) {
				mtx.val[i*mtx.n + j] = -1;
				mtx.val[j*mtx.n + i] = -1;
				mtx.val[j*mtx.n + j] += 1;
				mtx.val[i*mtx.n + i] += 1;
			}
		}
	}

	// console.log("Calculated " + ((mtx.n*(mtx.n+1))/2) + " distances in: " + (t1 - t0) + " milliseconds");
	
	// Calculate normal modes (use Lalolib)
	t0 = performance.now();
	E = eig(mtx, true);
	t1 = performance.now();
	console.log("Diagonalized matrix in: " + (t1 - t0) + " milliseconds");
	
	// Filter zero-value modes
	const small = 0.00001; // 1e-5
	const numZeroModes = sum(E.V.map(eval => (eval <= small) ? 1: 0));
	console.log(numZeroModes + " zero modes detected")

	sortedIndices = sort(E.V, false, true);
	eigenvalues = E.V.slice(numZeroModes, numModes + numZeroModes);
	const euT = transpose(E.U);
	eigenvectors = sortedIndices.slice(numZeroModes, numModes + numZeroModes).map(i => euT.row(i));

	// extendModes();
	console.log(eigenvalues);
	modeVariance = eigenvalues.map(ev => (1/ev));

	let selectedModes = [0, 1];
	const atomicFluctuations = calcSquareFluctuations(selectedModes);
	console.log(atomicFluctuations);

	return mtx
}

//
// Auxiliary functions
//
function sqEuclideanDistance(a, b) {
	return ((a[0] - b[0])*(a[0] - b[0]) + 
		    (a[1] - b[1])*(a[1] - b[1]) + 
		    (a[2] - b[2])*(a[2] - b[2]))
}

function calcSquareFluctuations(modeIndices) {

	if (modeIndices.length == 1) { // only one
		const m = modeIndices[0];
		return entrywisemul(pow(eigenvectors[m], 2), modeVariance[m]);
	}
	else {
		let squareFluctuations = zeros(numAtoms);
		for (let i = 0; i < modeIndices.length; i++) {
			const m = modeIndices[i];
			const modeFluct = entrywisemul(pow(eigenvectors[m], 2), modeVariance[m]);
			squareFluctuations = squareFluctuations.map((num, idx) => num + modeFluct[idx]);
		}
		return squareFluctuations;
	}
}

// Update coordinates
function updateCoords() {
	let structure = component.structure;

	const p = { what: { position: true } };
	const initialCoords = structure.getAtomData(p).position;
	const shuffledCoords = initialCoords.map(x => x + Math.round(Math.random() * 2.0));

	structure.updatePosition(shuffledCoords);
	component.updateRepresentations({ 'position': true });
}

function extendModes() {
	// Extends normal modes from Ca to all-atoms in the structure
	// Assumes CA order is not changed

	let ext_eigenval = [];
	let ext_eigenvec = [];

	const structureObject = component.structure;
	let atomCounter = 0;
	structureObject.eachResidue(function(res) {
		if (res.hasAtomWithName("CA")) {
			for (let i = 0; i < res.atomCount; i++) {
				ext_eigenval.push(eigenvectors[atomCounter]);
				ext_eigenvec.push(eigenvalues[atomCounter]);
			}
			atomCounter++;
		};
	});
	eigenvectors = ext_eigenvec;
	eigenvalues = ext_eigenval;
}

// Make lineplot of atomic fluctuations
// Make bar plot of eigenvectors
// Plot arrow vectors on structure (or) animate coordinates