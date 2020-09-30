import { to3ptModel } from "./cgmodels.js"
import { GNM } from './nma.js';
import { downloadMolecule } from './parser.js';
import { displayPDB } from './viewer.js';

// Bind function to key presses
let progressbox = document.getElementById("progress-text");
let textbox = document.getElementById("input-text");
textbox.addEventListener("keydown", function (e) {
    if (e.key === "Enter") {  //checks whether the pressed key is "Enter"
        // Hide box and show progress text
        textbox.className = "hidden";
        progressbox.className = "visible";
        // Run calculation
        run(e);
    }
});

export async function run(e) {

    progressbox.innerText = "running"
    let nmodes = 20;

    const pdbcode = e.target.value;

    const pdbdata = await downloadMolecule(pdbcode);

    const cgmodel = to3ptModel(pdbdata);

    GNM(cgmodel, nmodes);

    displayPDB(cgmodel);
}
