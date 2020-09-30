/**
 * Builds a connectivity (Kirchhoff) matrix from coordinates.
 * @param {Float32Array} coords - flat array of 3D coordinates.
 * @param {number} threshold - max distance to consider two atoms connected (in A).
 */
export function buildKirchhoff(coords, threshold=10.0) {

    const dt = threshold * threshold;
    const natoms = coords.length / 3;

    let t0 = performance.now()
    // const ksize = (natoms * natoms - natoms) / 2
    const ksize = natoms * natoms;
    const k = new Float32Array(ksize);

    // We use a flattened array for each RESIDUE
    // contact, but we iterate on a 3xN array (xyz/resid)
    for (let i = 0; i < coords.length - 3; i += 3) {

        let xi = coords[i];
        let yi = coords[i + 1];
        let zi = coords[i + 2];

        let si = i / 3;  // stride in k array

        for (let j = i+3; j < coords.length; j += 3) {

            let xj = coords[j];
            let yj = coords[j + 1];
            let zj = coords[j + 2];

            let dij = sqEuclideanDistance(xi, yi, zi, xj, yj, zj);
            if (dij <= dt) {

                let sj = j / 3;

                // off diagonal
                k[si * natoms + sj] = -1;
                k[sj * natoms + si] = -1;
                // diagonal terms
                k[si * natoms + si] += 1;
                k[sj * natoms + sj] += 1;
            }
        }
    }

    let t1 = performance.now()
    console.log('Built K matrix for ' + natoms + ' atoms in ' + (t1 - t0) + ' ms');
    return [natoms, k];
}

function sqEuclideanDistance(xi, yi, zi, xj, yj, zj) {
    return (
        (xi - xj) * (xi - xj) +
        (yi - yj) * (yi - yj) +
        (zi - zj) * (zi - zj)
    )
}

// Fluctuations from eigenvalues/vectors
