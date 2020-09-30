
const fullCode = new RegExp("^[0-9A-Za-z]{4}$");
const chainCode = new RegExp("([0-9A-Za-z]{4})\.([A-Z0-9a-z]+)$");

/**
 * Fetches a molecule from RCSB PDB.
 * @param {string} code - 4-letter PDB identifier
 */
export async function downloadMolecule(code) {
    console.log('Downloading pdb code: ' + code)

    // Do we have 1ctf or 1ctf.A
    let pdbcode;
    let chains;
    if (chainCode.test(code)) { // e.g. 1ctf.A
        let match = chainCode.exec(code);
        pdbcode = match[1];
        chains = match[2].split("");
    }
    else if (fullCode.test(code)) {
        pdbcode = code;
    }
    else {
        alert("Wrong input: pdb code must be e.g. 1ctf or 1ctf.A where A is a chain ID.")
        throw "wrong input: " + code;
    }

    let url = "https://files.rcsb.org/view/" + pdbcode + ".pdb";
    let response = await fetch(url);
    return filterChains(await response.text(), chains);
}

async function filterChains(data, chainlist) {

    const pdblines = data.split("\n");
    const chainset = new Set(chainlist);

    let lines = [];
    for (let i = 0, nlines = data.length; i < nlines; ++i) {

        let line = pdblines[i];
        if (line == null) {
            break  // EoF
        }
        if (line.startsWith("ATOM") && (chainset.size == 0 || chainset.has(line[21]))) {
            lines.push(line)
        }
    }
    return lines
}

/**
 * Extracts 3D coordinates from the raw PDB data.
 * @param {string} pdbdata - PDB file as text.
 */
// export function getAtomCoord(pdbdata) {

//     var t0 = performance.now()

//     let pdblines = pdbdata.split("\n");
//     let rawCoords = [];

//     for (let i = 0; i < pdblines.length; i += 1) {

//         if (pdblines[i].startsWith("ATOM")) {

//             let name = pdblines[i].slice(12, 16).replace(/\s/g, "");
//             if (name != "CA") {  // Only CA atoms
//                 continue
//             }

//             let x = Number(pdblines[i].slice(30, 38))
//             let y = Number(pdblines[i].slice(38, 46))
//             let z = Number(pdblines[i].slice(46, 54))

//             rawCoords.push(x);
//             rawCoords.push(y);
//             rawCoords.push(z);
//         }
//     }

//     let coords = Float32Array.from(rawCoords);

//     var t1 = performance.now()
//     console.log('Parsed ' + coords.length + ' coordinates in ' + (t1 - t0) + ' ms');
//     return coords;
// }
