import { constants } from "./constants.js";
import { getStructure } from './structure.js';
import { buildDistanceMatrix } from "./geometry.js";
// import { GNM } from './nma.js';
// import { displayPDB } from './viewer.js';


// Bind function to key presses
let progressContainer = document.getElementById("progress-text");
let inputContainer = document.getElementById("input-container");
let inputTextBox = document.getElementById("input-text");
inputTextBox.addEventListener("keydown", function (e) {
    if (e.key === "Enter") {  //checks whether the pressed key is "Enter"

        // Hide box and show progress text
        inputContainer.className = "hidden";
        progressContainer.className = "visible";

        // Run calculation
        run(e);
    }
});

// GL Viewer Things
constants.glviewer = $3Dmol.createViewer("gldiv", {});  // global scope

// "Main" function
export async function run(e) {

    const pdbcode = e.target.value;

    const structure = await getStructure(pdbcode);
    constants.distmtx = buildDistanceMatrix(structure);  // set globally _once_

    // displayPDB(cgmodel);
}
