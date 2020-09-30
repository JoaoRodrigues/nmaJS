// Mapping between side-chain atoms and 3pt COM beads.
// COM beads are placed on this particular atom or on the center of mass
// if there are multiple.
// Reference:
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1802011/
//
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

const centralBeads = new Set(["C4*", "CA"])  // for bonding

// Utility function to test if an object is empty.
function isEmptyObj(o) {
    return Object.keys(o).length === 0;
}

// Utility function to deep copy object
function copyObject(o) {
    return JSON.parse(JSON.stringify(o));
}

/**
 * Converts an all-atom structure to a 3pt model for 3DMol.js
 * @param {string} pdbdata - PDB file as text.
 */
export function to3ptModel(pdbdata) {

    const pdblines = pdbdata;

    let cgbeads = [];

    let lastResidue = null; // unique residue ID
    let resLines = [];  // holds each residue's lines.
    for (let i = 0, nlines = pdblines.length; i < nlines; ++i) {

        let line = pdblines[i]

        if (!line.startsWith("ATOM")) {
            continue  // we ignore HETATMs, sorry ligand people.
        }

        let residue = line.slice(17, 27);  // unique residue id.
        if (residue != lastResidue) {

            if (resLines.length) {
                let beads = residue2bead(resLines);
                cgbeads.push(...beads);

                resLines = [];
            }
            lastResidue = residue;
        }

        resLines.push(line)
    }

    // Catch Last
    if (resLines.length) {
        let beads = residue2bead(resLines);
        cgbeads.push(...beads);
    }

    // Assign serial numbers and bonds between consecutive central beads
    let lastCentral = {};
    let centralMap = {};  // for intra-residue bonds
    for (let i = 0, nbeads = cgbeads.length; i < nbeads; ++i) {
        let bead = cgbeads[i];
        bead.serial = i;

        if (centralBeads.has(bead.atom)) {
            if (bead.resi != lastCentral.resi && bead.chain == lastCentral.chain) {
                bead.bonds.push(lastCentral.serial);
            }
            lastCentral = bead;
            centralMap[bead.resi + '_' + bead.chain] = bead.serial
        }
    }

    // Now assign intra-residue bonds (mostly for viz)
    // CA to all in case of amino
    // C4* to all in case of nucleic
    for (let i = 0, nbeads = cgbeads.length; i < nbeads; ++i) {
        let bead = cgbeads[i];

        if (centralBeads.has(bead.atom)) {
            continue
        }

        let idx = centralMap[bead.resi + '_' + bead.chain];
        cgbeads[idx].bonds.push(i)
    }

    console.log('Created ' + cgbeads.length + ' beads');
    return cgbeads;
}

/**
 * Converts an all-atom residue into a CG bead.
 * @param {string} residue - residue data as PDB text format.
 */
function residue2bead(residue) {

    const resn = residue[0].slice(17, 20).replace(/\s/g, "");
    const resi = Number(residue[0].slice(22, 26).replace(/\s/g, ""));
    const chain = residue[0].slice(21, 22).replace(/\s/g, "");

    if (!(resn in scbeads)) {
        console.log(`Unsupported residue: ${resn}`)
        return [];
    }

    let beads = [];
    let comxx = 0;
    let comyy = 0;
    let comzz = 0;
    let nscat = 0;

    const scbeadmap = scbeads[resn];
    const defbead = {
        resn: resn,
        resi: resi,
        chain: chain,
        bonds: []
    }

    let centralbead = null;  // index of central bead in beads array.

    for (let i = 0, natoms = residue.length; i < natoms; ++i) {

        let name = residue[i].slice(12, 16).replace(/\s/g, "");
        let x = Number(residue[i].slice(30, 38))
        let y = Number(residue[i].slice(38, 46))
        let z = Number(residue[i].slice(46, 54))

        if (scbeadmap.has(name)) {
            comxx += x;
            comyy += y;
            comzz += z;
            nscat += 1;
        }
        else if (name == "CA") {
            let bead = copyObject(defbead);
            bead.atom = name
            bead.x = x
            bead.y = y
            bead.z = z
            bead.elem = "C"

            beads.push(bead)
            centralbead = beads.length - 1;
        }
        else if (name == "O") {
            let bead = copyObject(defbead);
            bead.atom = name
            bead.x = x
            bead.y = y
            bead.z = z
            bead.elem = "O"

            // beads.push(bead)
        }
        else if (name == "C4*" || name == "C4'") {
            let bead = copyObject(defbead);
            bead.atom = "C4*"
            bead.x = x
            bead.y = y
            bead.z = z
            bead.elem = "C"

            beads.push(bead)
            centralbead = beads.length - 1;
        }
        else if (name == "P") {
            let bead = copyObject(defbead);
            bead.atom = name
            bead.x = x
            bead.y = y
            bead.z = z
            bead.elem = "P"

            beads.push(bead)
        }
    }

    if (centralbead === null) {
        console.log(`Residue missing central bead: ${chain}:${resn}${resi}`)
        return [];
    }

    // Now process SC bead if necessary
    if (nscat > 0) {

        if (Number.isNaN(comxx) || Number.isNaN(comyy) || Number.isNaN(comzz)) {
            console.log(`COM coordinates undefined for residue: ${chain}:${resn}${resi}`)
            console.log(`comxx = ${comxx}, comyy = ${comyy}, comzz = ${comzz}, nscat = ${nscat}`)
            return [];
        }

        comxx /= nscat;
        comyy /= nscat;
        comzz /= nscat;

        let bead = copyObject(defbead);
        bead.atom = "COM"  // these are called CMA/CMN in the original paper
        bead.x = comxx
        bead.y = comyy
        bead.z = comzz
        bead.elem = "C"  // why not..

        // beads.push(bead)
    }

    // Special case for glycine: CA is also CMA
    if (resn == "GLY") {
        let cabead = beads[centralbead];
        let bead = copyObject(cabead);
        bead.atom = "COM"

        // beads.push(bead)
    }

    return beads;
}
