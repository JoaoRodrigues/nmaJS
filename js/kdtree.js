// Adapted from NGL.js
// Removed maxNodes parameter in nearest function

"use strict";
/**
 * @file Kdtree
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @private
 */
/**
 * Kdtree
 * @class
 * @author Alexander Rose <alexander.rose@weirdbyte.de>, 2016
 * @author Roman Bolzern <roman.bolzern@fhnw.ch>, 2013
 * @author I4DS http://www.fhnw.ch/i4ds, 2013
 * @license MIT License <http://www.opensource.org/licenses/mit-license.php>
 * @description
 * k-d Tree for typed arrays of 3d points (e.g. for Float32Array), in-place
 * provides fast nearest neighbour search
 *
 * Based on https://github.com/ubilabs/kd-tree-javascript by Ubilabs
 *
 * Further information (including mathematical properties)
 * http://en.wikipedia.org/wiki/Binary_tree
 * http://en.wikipedia.org/wiki/K-d_tree
 *
 * @example
 * points: [x, y, z, x, y, z, x, y, z, ...]
 * metric: function(a, b){
 *    return Math.pow(a[0]-b[0], 2) + Math.pow(a[1]-b[1], 2) + Math.pow(a[2]-b[2], 2);
 * }
 *
 * @param {Float32Array} points - points
 * @param {Function} metric - metric
 */
export class KDTree {
    constructor(points, metric) {
        this.points = points;
        this.metric = metric;
        this.maxDepth = 0;
        this.currentNode = 0;
        const n = points.length / 3;
        const indices = new Uint32Array(n);
        for (let i = 0; i < n; ++i) {
            indices[i] = i;
        }
        this.indices = indices;
        this.nodes = new Int32Array(n * 4);
        this.rootIndex = this.buildTree(0, -1, 0, n);
    }
    buildTree(depth, parent, arrBegin, arrEnd) {
        if (depth > this.maxDepth)
            this.maxDepth = depth;
        const plength = arrEnd - arrBegin;
        if (plength === 0) {
            return -1;
        }
        const nodeIndex = this.currentNode * 4;
        const nodes = this.nodes;
        this.currentNode += 1;
        if (plength === 1) {
            nodes[nodeIndex] = arrBegin;
            nodes[nodeIndex + 1] = -1;
            nodes[nodeIndex + 2] = -1;
            nodes[nodeIndex + 3] = parent;
            return nodeIndex;
        }
        // if(plength <= 32){
        //   return nodeIndex;
        // }
        const indices = this.indices;
        const points = this.points;
        const arrMedian = arrBegin + Math.floor(plength / 2);
        const currentDim = depth % 3;
        // inlined quickselect function
        let j, tmp, pivotIndex, pivotValue, storeIndex;
        let left = arrBegin;
        let right = arrEnd - 1;
        while (right > left) {
            pivotIndex = (left + right) >> 1;
            pivotValue = points[indices[pivotIndex] * 3 + currentDim];
            // swap( pivotIndex, right );
            tmp = indices[pivotIndex];
            indices[pivotIndex] = indices[right];
            indices[right] = tmp;
            storeIndex = left;
            for (j = left; j < right; ++j) {
                if (points[indices[j] * 3 + currentDim] < pivotValue) {
                    // swap( storeIndex, j );
                    tmp = indices[storeIndex];
                    indices[storeIndex] = indices[j];
                    indices[j] = tmp;
                    ++storeIndex;
                }
            }
            // swap( right, storeIndex );
            tmp = indices[right];
            indices[right] = indices[storeIndex];
            indices[storeIndex] = tmp;
            pivotIndex = storeIndex;
            if (arrMedian === pivotIndex) {
                break;
            }
            else if (arrMedian < pivotIndex) {
                right = pivotIndex - 1;
            }
            else {
                left = pivotIndex + 1;
            }
        }
        nodes[nodeIndex] = arrMedian;
        nodes[nodeIndex + 1] = this.buildTree(depth + 1, nodeIndex, arrBegin, arrMedian);
        nodes[nodeIndex + 2] = this.buildTree(depth + 1, nodeIndex, arrMedian + 1, arrEnd);
        nodes[nodeIndex + 3] = parent;
        return nodeIndex;
    }
    getNodeDepth(nodeIndex) {
        const parentIndex = this.nodes[nodeIndex + 3];
        return (parentIndex === -1) ? 0 : this.getNodeDepth(parentIndex) + 1;
    }
    // TODO
    // function getNodePos (node) {}
    /**
     * find nearest points
     * @param {Array} point - array of size 3
     * @param {Float} maxDistance - maximum distance of point to result nodes
     * @return {Array} array of point, distance pairs
     */
    nearest(point, maxDistance) {
        const bestNodes = new BinaryHeap(e => -e[1]);
        const nodes = this.nodes;
        const points = this.points;
        const indices = this.indices;
        const nearestSearch = (nodeIndex) => {
            let bestChild, otherChild;
            const dimension = this.getNodeDepth(nodeIndex) % 3;
            const pointIndex = indices[nodes[nodeIndex]] * 3;
            const ownPoint = [
                points[pointIndex + 0],
                points[pointIndex + 1],
                points[pointIndex + 2]
            ];
            const ownDistance = this.metric(point, ownPoint);
            function saveNode(nodeIndex, distance) {
                bestNodes.push([nodeIndex, distance]);
            }
            const leftIndex = nodes[nodeIndex + 1];
            const rightIndex = nodes[nodeIndex + 2];
            // if it's a leaf
            if (rightIndex === -1 && leftIndex === -1) {
                if (ownDistance <= maxDistance) {
                    saveNode(nodeIndex, ownDistance);
                }
                return;
            }
            if (rightIndex === -1) {
                bestChild = leftIndex;
            }
            else if (leftIndex === -1) {
                bestChild = rightIndex;
            }
            else {
                if (point[dimension] <= points[pointIndex + dimension]) {
                    bestChild = leftIndex;
                }
                else {
                    bestChild = rightIndex;
                }
            }
            // recursive search
            nearestSearch(bestChild);
            if (ownDistance <= maxDistance) {
                saveNode(nodeIndex, ownDistance);
            }
        };
        nearestSearch(this.rootIndex);
        const result = [];
        for (let i = 0, il = bestNodes.size(); i < il; i += 1) {
            result.push(bestNodes.content[i]);
        }
        return result;
    }
    verify(nodeIndex, depth = 0) {
        let count = 1;
        if (nodeIndex === undefined) {
            nodeIndex = this.rootIndex;
        }
        if (nodeIndex === -1) {
            throw new Error('node is null');
        }
        const dim = depth % 3;
        const nodes = this.nodes;
        const points = this.points;
        const indices = this.indices;
        const leftIndex = nodes[nodeIndex + 1];
        const rightIndex = nodes[nodeIndex + 2];
        if (leftIndex !== -1) {
            if (points[indices[nodes[leftIndex]] * 3 + dim] >
                points[indices[nodes[nodeIndex]] * 3 + dim]) {
                throw new Error('left child is > parent!');
            }
            count += this.verify(leftIndex, depth + 1);
        }
        if (rightIndex !== -1) {
            if (points[indices[nodes[rightIndex]] * 3 + dim] <
                points[indices[nodes[nodeIndex]] * 3 + dim]) {
                throw new Error('right child is < parent!');
            }
            count += this.verify(rightIndex, depth + 1);
        }
        return count;
    }
}
/**
 * @file Binary Heap
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @private
 */
/**
 * Binary heap implementation
 * @class
 * @author http://eloquentjavascript.net/appendix2.htm
 * @param {Function} scoreFunction - the heap scoring function
 */
class BinaryHeap {
    constructor(scoreFunction) {
        this.scoreFunction = scoreFunction;
        this.content = [];
        this.scoreFunction = scoreFunction;
    }
    push(element) {
        // Add the new element to the end of the array.
        this.content.push(element);
        // Allow it to bubble up.
        this.bubbleUp(this.content.length - 1);
    }
    pop() {
        // Store the first element so we can return it later.
        const result = this.content[0];
        // Get the element at the end of the array.
        const end = this.content.pop();
        // If there are any elements left, put the end element at the
        // start, and let it sink down.
        if (end && this.content.length > 0) {
            this.content[0] = end;
            this.sinkDown(0);
        }
        return result;
    }
    peek() {
        return this.content[0];
    }
    remove(element) {
        const len = this.content.length;
        // To remove a value, we must search through the array to find it.
        for (let i = 0; i < len; i++) {
            if (this.content[i] === element) {
                // When it is found, the process seen in 'pop' is repeated
                // to fill up the hole.
                const end = this.content.pop();
                if (end && i !== len - 1) {
                    this.content[i] = end;
                    if (this.scoreFunction(end) < this.scoreFunction(element)) {
                        this.bubbleUp(i);
                    }
                    else {
                        this.sinkDown(i);
                    }
                }
                return;
            }
        }
        throw new Error('Node not found.');
    }
    size() {
        return this.content.length;
    }
    bubbleUp(n) {
        // Fetch the element that has to be moved.
        const element = this.content[n];
        // When at 0, an element can not go up any further.
        while (n > 0) {
            // Compute the parent element's index, and fetch it.
            const parentN = Math.floor((n + 1) / 2) - 1;
            const parent = this.content[parentN];
            // Swap the elements if the parent is greater.
            if (this.scoreFunction(element) < this.scoreFunction(parent)) {
                this.content[parentN] = element;
                this.content[n] = parent;
                // Update 'n' to continue at the new position.
                n = parentN;
            }
            else {
                // Found a parent that is less, no need to move it further.
                break;
            }
        }
    }
    sinkDown(n) {
        // Look up the target element and its score.
        const length = this.content.length;
        const element = this.content[n];
        const elemScore = this.scoreFunction(element);
        let child1Score = 0;
        let child2Score = 0;
        while (true) {
            // Compute the indices of the child elements.
            const child2N = (n + 1) * 2;
            const child1N = child2N - 1;
            // This is used to store the new position of the element, if any.
            let swap = null;
            // If the first child exists (is inside the array)...
            if (child1N < length) {
                // Look it up and compute its score.
                const child1 = this.content[child1N];
                child1Score = this.scoreFunction(child1);
                // If the score is less than our element's, we need to swap.
                if (child1Score < elemScore)
                    swap = child1N;
            }
            // Do the same checks for the other child.
            if (child2N < length) {
                const child2 = this.content[child2N];
                child2Score = this.scoreFunction(child2);
                if (child2Score < (swap === null ? elemScore : child1Score))
                    swap = child2N;
            }
            // If the element needs to be moved, swap it, and continue.
            if (swap !== null) {
                this.content[n] = this.content[swap];
                this.content[swap] = element;
                n = swap;
            }
            else {
                // Otherwise, we are done.
                break;
            }
        }
    }
}
