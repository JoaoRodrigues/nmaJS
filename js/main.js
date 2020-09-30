import { to3ptModel } from "./cgmodels.js"
import { GNM } from './nma.js';
import { downloadMolecule } from './parser.js';
import { displayPDB } from './viewer.js';

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

export async function run(e) {

    let nmodes = 20;

    const pdbcode = e.target.value;

    const pdbdata = await downloadMolecule(pdbcode);

    const cgmodel = to3ptModel(pdbdata);

    GNM(cgmodel, nmodes);

    displayPDB(cgmodel);
}
