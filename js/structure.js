import { constants } from "./constants.js";

const fullCode = new RegExp("^[0-9A-Za-z]{4}$");
const chainCode = new RegExp("([0-9A-Za-z]{4})\.([A-Z0-9a-z]+)$");

/**
 * Fetches and parses a molecule from RCSB PDB.
 * @param {string} code - 4-letter identifier (e.g. 1ctf) with
 *  optional chain IDs (e.g. 1ctf.A, 1brs.AD).
 */
export async function getStructure(code) {
    console.log('Downloading pdb code: ' + code)

    // Do we have 4-letter code or 4+chain?
    let pdbcode;
    let chains;
    if (chainCode.test(code)) { // e.g. 1ctf.A
        let match = chainCode.exec(code);
        pdbcode = match[1];
        chains = match[2].split("");
    }
    else if (fullCode.test(code)) {  // e.g. 1ctf
        pdbcode = code;
    }
    else {
        alert("Wrong input: pdb code must be e.g. 1xyz or 1xyz.A where A is the chain ID.")
        throw "wrong input: " + code;
    }

    let url = "https://mmtf.rcsb.org/v1.0/full/" + pdbcode + ".mmtf.gz";
    let response = await $3Dmol.getbin(url);
    return parseStructure(await response, chains);
}

/**
 * Uses 3DMol to parse the structure into an object.
 * @param {string} moldata - molecule data as a base64 string.
 * @param {array} chains - chain identifiers used to filter the structure.
 */
function parseStructure(moldata, chains) {

    let m = constants.glviewer.addModel();
    m.addMolData(moldata, "mmtf");

    filterChains(m, chains)
    filterHETATM(m)

    createCOM(m)

    return m  // return molecule
}

function filterChains(model, chainlist) {

    const chainset = new Set(chainlist);
    if (chainset.size == 0) { return }  // no selection

    let atoms = model.selectedAtoms({});
    let badatoms = [];
    for (let i = 0, natoms = atoms.length; i < natoms; ++i) {
        let atom = atoms[i];
        if (!chainset.has(atom.chain)) {
            badatoms.push(atom);
        }
    }

    model.removeAtoms(badatoms);
    console.log('Removed ' + badatoms.length + ' out of ' + atoms.length + ' atoms');
}

function filterHETATM(model) {

    let atoms = model.selectedAtoms({});
    let badatoms = [];
    for (let i = 0, natoms = atoms.length; i < natoms; ++i) {
        let atom = atoms[i];
        if (atom.hetflag) {
            badatoms.push(atom);
        }
    }

    model.removeAtoms(badatoms);
    console.log('Removed ' + badatoms.length + ' hetatms');
}

// Code to generate Center of Mass beads

/**
** Mapping between side-chain atoms and 3pt COM beads.
** COM beads are placed on this particular atom or on the center of mass
** if there are multiple.
** Reference:
** https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1802011/
**/
const scbeads = {
    // amino
    LEU: new Set(["CG"]),
    ALA: new Set(["CB"]),
    GLY: new Set([]),
    VAL: new Set(["CB"]),
    GLU: new Set(["CD"]),
    ASP: new Set(["CG"]),
    THR: new Set(["CB"]),
    LYS: new Set(["CD"]),
    ILE: new Set(["CG1"]),
    ARG: new Set(["NE"]),
    ASN: new Set(["CG"]),
    PRO: new Set(["CG"]),
    PHE: new Set(["CG"]),
    GLN: new Set(["CD"]),
    SER: new Set(["CB", "OG"]),
    HIS: new Set(["CG"]),
    MET: new Set(["CG"]),
    TYR: new Set(["CD1", "CD2"]),
    TRP: new Set(["CD2"]),
    CYS: new Set(["CB", "SG"]),

    // nucleic
    A: new Set(["C2", "C4", "C6"]),
    T: new Set(["C2", "C4", "C6"]),
    C: new Set(["C2", "C4", "C6"]),
    G: new Set(["C2", "C4", "C6"]),
    U: new Set(["C2", "C4", "C6"]),
    DA: new Set(["C2", "C4", "C6"]),
    DT: new Set(["C2", "C4", "C6"]),
    DC: new Set(["C2", "C4", "C6"]),
    DG: new Set(["C2", "C4", "C6"]),
}

/**
 * Creates center-of-mass beads for each residue.
 */
function createCOM(model) {

    let lastResidue = null;
    let resAtoms = [];  // to hold residue atoms.
    let beads = []; // to hold COM beds

    let atoms = model.selectedAtoms({});
    for (let i = 0, natoms = atoms.length; i < natoms; ++i) {
        let atom = atoms[i];

        // Residue unique id string
        let residue = atom.chain + atom.resi + atom.resn + atom.icode
        if (!!lastResidue && residue != lastResidue) {
            let bead = residue2bead(resAtoms);
            beads.push(...bead);
            resAtoms = [];  // clear
        }

        resAtoms.push(atom);
        lastResidue = residue;
    }

    // Catch last atom
    if (resAtoms.length) {
        let bead = residue2bead(resAtoms);
        beads.push(...bead);
        resAtoms = [];
    }

    console.log("Created " + beads.length + " center-of-mass beads.")
    model.addAtoms(beads);
}

/**
 * Converts an all-atom residue into a CG bead.
 * @param {array} atomlist - atom objects from 3DMol.js
 */
function residue2bead(atomlist) {

    const centralAtomSet = constants.centralAtoms;

    const resn = atomlist[0].resn;
    const resi = atomlist[0].resi;
    const chain = atomlist[0].chain;

    if (!(resn in scbeads)) {
        console.log("Unsupported residue for COM bead creation: " + resn)
        return [];
    }

    let centralatom = null;
    const scbeadmap = scbeads[resn];
    const scbead = {
        atom: "COM",
        resn: resn,
        resi: resi,
        chain: chain,
        elem: "C",  // random really..
        bonds: []
    }

    let comxx = 0;
    let comyy = 0;
    let comzz = 0;
    let nscat = 0;

    for (let i = 0, natoms = atomlist.length; i < natoms; ++i) {

        let atom = atomlist[i];

        let name = atom.atom;

        if (scbeadmap.has(name)) {
            comxx += atom.x;
            comyy += atom.y;
            comzz += atom.z;
            nscat += 1;
        }
        else if (centralAtomSet.has(name)) {
            centralatom = atom;
        }
    }

    if (centralatom === null) {
        console.log(`Residue missing central atom: ${chain}:${resn}${resi}`)
        return [];
    }

    // Special case for glycine: CA is also CMA
    if (resn == "GLY") {
        comxx += centralatom.x;
        comyy += centralatom.y;
        comzz += centralatom.z;
        nscat += 1;
    }

    // Now process SC bead if necessary
    if (nscat == 0) {
        console.log(`Residue ${chain}:${resn}${resi} does not have side-chain atoms.`)
        return [];
    }

    if (Number.isNaN(comxx) || Number.isNaN(comyy) || Number.isNaN(comzz)) {
        console.log(`COM coordinates undefined for residue: ${chain}:${resn}${resi}`)
        console.log(`comxx = ${comxx}, comyy = ${comyy}, comzz = ${comzz}, nscat = ${nscat}`)
        return [];
    }

    scbead.x = comxx / nscat;
    scbead.y = comyy / nscat;
    scbead.z = comzz / nscat;

    // Add bond to central bead
    scbead.bonds.push(centralatom.serial);
    // console.log(`Created COM bead for residue: ${chain}:${resn}${resi}`)

    return [scbead];  // we return as a list to be coherent with the empty returns
}
