import { constants } from "./constants.js";
import { KDTree } from "./kdtree.js"

/**
 * Builds a distance matrix from coordinates.
 * @param {object} structure - 3DMol.js GLModel object.
 * @param {number} threshold - max distance to consider (in A).
 */
export function buildDistanceMatrix(structure, threshold=5.0) {

    const dt = threshold * threshold;
    const centralAtomSet = constants.centralAtoms;
    let dim = 0; // number of central atoms (dimensions of dmtx)

    const atoms = structure.selectedAtoms({});

    // Extract coordinates and map atoms to residues.
    let atom2residue = {};
    let coords = [];
    for (let i = 0; i < atoms.length; ++i) {
        let atom = atoms[i];
        if (atom.atom == "COM") {
            continue // ignore centers of mass
        }

        coords.push(atom.x)
        coords.push(atom.y)
        coords.push(atom.z)

        let residue = atom.chain + atom.resi + atom.resn + atom.icode;
        atom2residue[i] = residue;

        if (centralAtomSet.has(atom.atom)) {
            dim++;
        }
    }
    coords = Float32Array.from(coords);
    const natoms = Math.floor(coords.length / 3);
    console.log(`Found ${dim} central atoms out of ${natoms} total atoms`);

    // Build KDTree
    let t0 = performance.now()
    const metric = sqEuclideanDistance;
    const kdt = new KDTree(coords, metric);
    const kdtIdx = kdt.indices;
    const kdtNodes = kdt.nodes;

    let t1 = performance.now()
    console.log('Built kdtree for ' + natoms + ' atoms in ' + (t1 - t0) + ' ms');

    // Build distance matrix.
    // We take the shortest distance between any two atoms of two residues.

    t0 = performance.now()
    for (let i = 0; i < natoms; ++i) {
        let atom = atoms[i];
        let residue = atom2residue[i];

        let neighbors = {};  // neighbor id -> shortest distance
        let result = kdt.nearest([atom.x, atom.y, atom.z], dt);
        for (let n = 0; n < result.length; ++n) {
            let nodeIdx = result[n][0];
            let nDist = result[n][1];
            let neighborIdx = kdtIdx[kdtNodes[nodeIdx]];
            let neighborResidue = atom2residue[neighborIdx];

            // Ignore same residue neighbors.
            if (neighborResidue == residue) {
                continue
            }

            if (neighborResidue in neighbors) {
                if (neighbors[neighborResidue] > nDist) {
                    continue
                }
            }
            neighbors[neighborResidue] = nDist;
            console.log(residue, atom.atom, neighborResidue, atoms[neighborIdx].atom, Math.sqrt(nDist))
        }
        // console.log(`Atom ${atom.atom} of residue ${residue} has ${neighbors.length} neighbors within ${threshold}A`)
        // console.log(atom, neighbors)
    }
    t1 = performance.now()
    console.log('Search completed in ' + (t1 - t0) + ' ms');

    // Populate distance matrix
    const dmtx = new Float32Array(dim * dim);
}

const sqEuclideanDistance = function(a, b) {
    const dx = a[0] - b[0];
    const dy = a[1] - b[1];
    const dz = a[2] - b[2];

    return dx * dx + dy * dy + dz * dz;
}















// /**
//  * Builds a connectivity (Kirchhoff) matrix from coordinates.
//  * @param {Float32Array} coords - flat array of 3D coordinates.
//  * @param {number} threshold - max distance to consider two atoms connected (in A).
//  */
// export function buildKirchhoff(coords, threshold=10.0) {

//     const dt = threshold * threshold;
//     const natoms = coords.length / 3;

//     let t0 = performance.now()
//     // const ksize = (natoms * natoms - natoms) / 2
//     const ksize = natoms * natoms;
//     const k = new Float32Array(ksize);

//     // We use a flattened array for each RESIDUE
//     // contact, but we iterate on a 3xN array (xyz/resid)
//     for (let i = 0; i < coords.length - 3; i += 3) {

//         let xi = coords[i];
//         let yi = coords[i + 1];
//         let zi = coords[i + 2];

//         let si = i / 3;  // stride in k array

//         for (let j = i+3; j < coords.length; j += 3) {

//             let xj = coords[j];
//             let yj = coords[j + 1];
//             let zj = coords[j + 2];

//             let dij = sqEuclideanDistance(xi, yi, zi, xj, yj, zj);
//             if (dij <= dt) {

//                 let sj = j / 3;

//                 // off diagonal
//                 k[si * natoms + sj] = -1;
//                 k[sj * natoms + si] = -1;
//                 // diagonal terms
//                 k[si * natoms + si] += 1;
//                 k[sj * natoms + sj] += 1;
//             }
//         }
//     }

//     let t1 = performance.now()
//     console.log('Built K matrix for ' + natoms + ' atoms in ' + (t1 - t0) + ' ms');
//     return [natoms, k];
// }

// function sqEuclideanDistance(xi, yi, zi, xj, yj, zj) {
//     return (
//         (xi - xj) * (xi - xj) +
//         (yi - yj) * (yi - yj) +
//         (zi - zj) * (zi - zj)
//     )
// }

// /**
//  * Builds a connectivity (Kirchhoff) matrix from coordinates.
//  * Defines residues as neighbors if they have at least two
//  * heavy atoms within <threshold> of each other.
//  * @param {Float32Array} coords - flat array of 3D coordinates.
//  * @param {number} threshold - max distance to consider two atoms connected (in A).
//  */
// export function buildKirchhoff_v2(coords, threshold=5.0) {

//     const dt = threshold * threshold;
//     const natoms = coords.length / 3;

//     let t0 = performance.now()
//     // We use a flattened array for each RESIDUE
//     // contact, but we iterate on a 3xN array (xyz/resid)
//     for (let i = 0; i < coords.length - 3; i += 3) {

//         let pt = {"x": coords[i], "y": coords[i + 1], "z": coords[i + 2]};
//         let si = i / 3;  // stride in k array
//         let neighbors = tree.nearest(pt, 100, dt);
//         console.log(neighbors)
//     }

//     let t1 = performance.now()
//     console.log('Built K matrix for ' + natoms + ' atoms in ' + (t1 - t0) + ' ms');
//     // return [natoms, k];
// }

// /**
//  * Builds a connectivity (Kirchhoff) matrix from coordinates.
//  * Uses kdtree
//  * @param {Float32Array} coords - flat array of 3D coordinates.
//  * @param {number} threshold - max distance to consider two atoms connected (in A).
//  */
// export function buildKirchhoff_kdtree(coords, threshold=10.0) {

//     const dt = threshold * threshold;
//     const natoms = coords.length / 3;

//     // Convert points to array of objects
//     let points = [];
//     for (let i = 0; i < coords.length; i += 3) {
//         points.push(
//             {"x": coords[i], "y": coords[i + 1], "z": coords[i + 2]}
//         )
//     }
//     console.log('Parsed ' + points.length + ' points');

//     // Build kdTree
//     let t0 = performance.now()
//     const tree = new kdTree(points, PointSqDistance, ["x", "y", "z"]);
//     let t1 = performance.now()
//     console.log('Built kdtree for ' + natoms + ' atoms in ' + (t1 - t0) + ' ms');

//     t0 = performance.now()
//     // We use a flattened array for each RESIDUE
//     // contact, but we iterate on a 3xN array (xyz/resid)
//     for (let i = 0; i < coords.length - 3; i += 3) {

//         let pt = {"x": coords[i], "y": coords[i + 1], "z": coords[i + 2]};
//         let si = i / 3;  // stride in k array
//         let neighbors = tree.nearest(pt, 100, dt);
//         console.log(neighbors)
//     }

//     t1 = performance.now()
//     console.log('Built K matrix for ' + natoms + ' atoms in ' + (t1 - t0) + ' ms');
//     // return [natoms, k];
// }

// const PointSqDistance = function(a, b) {
//     return (
//         (a.x - b.x) * (a.x - b.x) +
//         (a.y - b.y) * (a.y - b.y) +
//         (a.z - b.z) * (a.z - b.z)
//     )
// }
