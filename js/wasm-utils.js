// WASM utility functions

// export function arrayToHeap(typedArray) {
//     console.log('Allocating: ' + typedArray.length + ' * ' + typedArray.BYTES_PER_ELEMENT + ' bytes');

//     let databytes = typedArray.length * typedArray.BYTES_PER_ELEMENT;
//     let dataptr = Module._malloc(databytes);

//     let heapbytes = new Uint8Array(Module.HEAPU8.buffer, dataptr, databytes);
//     heapbytes.set(new Uint8Array(typedArray.buffer));

//     return heapbytes;
// }

/**
 * Helper function to share arrays between JS and WASM.
 * @param {Float32Array} F32Array - array of floats to pass to WASM code.
 */
export function F32toHeap(F32Array) {
    // console.log('Allocating: ' + F32Array.length + ' * ' + F32Array.BYTES_PER_ELEMENT + ' bytes');

    let buffer = Module._malloc(F32Array.length * F32Array.BYTES_PER_ELEMENT);
    Module.HEAPF32.set(F32Array, buffer >> 2);  // copy the array

    return buffer
}
