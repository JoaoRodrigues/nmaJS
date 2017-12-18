// Main JavaScript Script

// Global vars
let stage;
let component;

// UI Global Variables
let atomPicked = null;
let bondPicked = null;
let bondsDrawn = new Set();

// Regexes
const fullStructure = new RegExp("[0-9A-Za-z]{4}$");
const partialStructure = new RegExp("([0-9A-Za-z]{4})\.([A-Z0-9a-z]{1})$");

// RGB Arrays for colors
const newBondColor = "rgb(255, 155, 0)" // Light orange
const selectedHaloColor = "rgb(255, 255, 150)" // Yellowish
const highlightedBondColor = "rgb(200, 0, 0)"; // Red

async function launchNMA(e) {

    // Track keyboard events until 'Enter' is pressed (code=13)
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

    // Get structure from CoordinateServer
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
    
    // Load Structure
    await updateProgress("Fetching your favorite protein");
    try {
        component = await stage.loadFile(url, { ext: "cif" });  
    }
    catch (err) {
        await updateProgress("Are you sure *this* is your favorite protein: " + structureCode);
        return false;
    }

    await updateProgress("Making it look pretty");
    await renderMolecule();  // Render molecule as hyperballs. Return structure object

    // Run function from external module
    await buildKirchhoff(component.structure);
    await NMA(mtx);
    //

    await colorByFluctuation(sqfluctuations);

    // Finally, fade out load box
    await updateProgress("Ta-dah! Have fun!", 1000);

    document.getElementById("loadbox").className = "hidden";
    document.getElementById("progress_div").className = "invisible"

    // Update stats div
    let structureNumbers = "PDB ID: " + structureCode + " | No. Atoms: " + component.structure.atomCount + " | No. Network Nodes: " + numAtoms;
    document.getElementById("stats").innerHTML = structureNumbers;

    return true;

}

function renderMolecule() {

    component.addRepresentation('hyperball', { sele: 'protein', 
                                               color: 'white',
                                               quality: 'high' }); 

    let pa = component.structure.getView(new NGL.Selection('.CA')).getPrincipalAxes();
    stage.animationControls.rotate(pa.getRotationQuaternion(), 0);
    stage.autoView();

    stage.mouseControls.add("clickPick-left", pickElement);
    stage.mouseControls.add("hoverPick", hoverComponent);
}

// Updates messages under input text with a certain delay (in miliseconds)
function updateProgress(message, delay) {

    let el = document.getElementById('progress_div');
    el.innerHTML = message;

    // Hackish way to allow DOM/UI updates
    return new Promise((resolve, reject) => {
        setTimeout(() => resolve(true), delay);
    });
}

// UI functions
// Most important functions of *this* universe

// For simplicity, let's use a single-button UI
// One click on an atom creates a selection (and an accompanying halo)
// A second click on another atom draws a bond and clears the selection.
// Clicking a halo undoes the selection.
// Clicking a bond twice removes it (first time highlights red to convey selection)
function pickElement(stage, pickingProxy) {

    // Click on Halo
    if (pickingProxy && pickingProxy.sphere) {
        if (pickingProxy.sphere.shape.name == "SelectedAtomHalo") {
            // Clear selection and remove halo
            atomPicked = null;
            destroyHalo();
        }
    }

    // Click on Atom
    if (pickingProxy && (pickingProxy.atom || pickingProxy.closestBondAtom)) {
        let atom = pickingProxy.atom || pickingProxy.closestBondAtom;

        if (atomPicked == null) { // first atom picked
            atomPicked = atom;
            drawHalo(atom, selectedHaloColor);
        }
        else {

            let other = atomPicked;
            
            // This should not be possible because the halo sphere overlaps with the 
            // atom representation. But just for sanity.
            if (other != atom) {
                // If bond does not exist, draw bond and calculate changes
                let bondName = "bond_" + atom.index + "_" + other.index;
                if (!bondsDrawn.has(bondName)) {

                    let iA = residueToIdx[other.resno];
                    let iB = residueToIdx[atom.resno];                
                    editKirchhoff(iA, iB, 10);
                    NMA(mtx);
                    colorByFluctuation(sqfluctuations);

                    drawBond(atom, other, newBondColor);
                }
            }

            // Clear selection and remove halo
            atomPicked = null;
            destroyHalo();
        } 
    };

    // Click on Bond
    if (pickingProxy && pickingProxy.cylinder) {

        let bond = pickingProxy.cylinder;
        let bondName = pickingProxy.cylinder.shape.name;

        // Get bonded atoms from bond name
        let tokens = bondName.split("_");
        let atomA = component.structure.getAtomProxy(tokens[1]);
        let atomB = component.structure.getAtomProxy(tokens[2]);

        if (bondPicked == null) { // Highlight bond
            bondPicked = bond;
            // Easier to redraw bond
            destroyBond(atomA, atomB);
            drawBond(atomA, atomB, highlightedBondColor);
        }
        else if (bondPicked.shape.name == bond.shape.name) { // Remove bond
            // Recalculate NMA
            let iA = residueToIdx[atomA.resno];
            let iB = residueToIdx[atomB.resno];
            editKirchhoff(iA, iB, -1); // negative factor restores original matrix element
            NMA(mtx);
            colorByFluctuation(sqfluctuations);
            // Destroy Bond
            destroyBond(atomA, atomB);

            bondPicked = null;
        }
        else { // Another bond was selected. Change selection.
            let otherBond = bondPicked;
            let otherBondName = otherBond.shape.name;
            let otherAtomA = component.structure.getAtomProxy(otherBondName.split("_")[1]);
            let otherAtomB = component.structure.getAtomProxy(otherBondName.split("_")[2]);

            destroyBond(otherAtomA, otherAtomB);
            drawBond(otherAtomA, otherAtomB, newBondColor);

            destroyBond(atomA, atomB);
            drawBond(atomA, atomB, highlightedBondColor);

            bondPicked = bond;
        }
    }
};

function drawHalo(atom, color = "red") {
    let haloColor = new NGL.Color(color);
    let shape = new NGL.Shape("SelectedAtomHalo");
    shape.addSphere([atom.x, atom.y, atom.z], haloColor, 1.0);
    let shapeComponent = stage.addComponentFromObject(shape);
    shapeComponent.addRepresentation("buffer");    
}

// Since we only make *one* halo, no need to get fancy
function destroyHalo() {
    let halo = stage.getComponentsByName("SelectedAtomHalo").first; // returns a collection with one element
    stage.removeComponent(halo);
}

function drawBond(atomA, atomB, color = "yellow") {

    let bondColor = new NGL.Color(color);

    let bondName = "bond_" + atomA.index + "_" + atomB.index;
    bondsDrawn.add(bondName);
    let rev_bondName = "bond_" + atomB.index + "_" + atomA.index;
    bondsDrawn.add(rev_bondName);

    let shape = new NGL.Shape(bondName, { openEnded: true }); // not capped
    shape.addCylinder([atomA.x, atomA.y, atomA.z], [atomB.x, atomB.y, atomB.z], bondColor, 0.2);

    let shapeComponent = stage.addComponentFromObject(shape);
    shapeComponent.addRepresentation("buffer");
}

function destroyBond(atomA, atomB) {

    let bondName = "bond_" + atomA.index + "_" + atomB.index;
    bondsDrawn.delete(bondName);
    let rev_bondName = "bond_" + atomB.index + "_" + atomA.index;
    bondsDrawn.delete(rev_bondName);

    let bond = stage.getComponentsByName(bondName).first; // returns a collection with one element
    stage.removeComponent(bond);
}

// Hovering
function hoverComponent(stage, pickingProxy) {

    let tooltip = document.getElementById("ngl-tooltip");

    if (pickingProxy && (pickingProxy.atom || pickingProxy.bond)){
    
        let atom = pickingProxy.atom || pickingProxy.closestBondAtom;
        let cp = pickingProxy.canvasPosition;

        let atomFluctuation = sqfluctuations[atom.index];
        tooltip.innerText = atom.resname + atom.resno + "." + atom.atomname + " = " + atomFluctuation.toFixed(3);
    
        tooltip.style.bottom = cp.y + 3 + "px";
        tooltip.style.left = cp.x + 3 + "px";
        tooltip.style.display = "block";
    }
    else {
        tooltip.style.display = "none";
    }
};

// Coloring Function
function colorByFluctuation(fluctarray) {
    // Use D3.js scale chromatic package
    let colorScale = d3.interpolateRdBu;

    // Normalize fluctuations
    // const maxSqF = Math.max.apply(null, fluctarray);
    // const minSqF = Math.min.apply(null, fluctarray);
    const minSqF = 0.0;
    const rangeSqF = (initialMaximumFluctuation - minSqF);
    let normalizedArray = fluctarray.map(v => ((v - minSqF) / rangeSqF));

    let schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
      this.atomColor = function (atom) {
        if (atom.hetero == true) {
            return 0x000000 // white
        } else {
            let rgbString = colorScale(1 - normalizedArray[atom.index]); // Reverse color scale the lazy way
            let col3 = d3.rgb(rgbString); // convert from string to rgb tuple
            let hexCode = col3 ? (col3.r << 16) + (col3.g << 8) + col3.b : 255; // rgb to hex (credit to martingraham)
            return hexCode;
        }
      }
    })

    component.reprList[0].setColor(schemeId);
}