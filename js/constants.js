// Module to define constants and abuse them
// to define mutable global variables between modules.

export const constants = {
    // 3DMol.js parameters
    "glviewer": null,

    // Structure/Geometry
    "centralAtoms": new Set(["C4*", "C4'", "CA"]),

    // GNM parameters
    "distmtx": null,
    "nmodes": 20,
}
