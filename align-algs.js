// Get user input and run the requested algorithm
seqAlignForm.onsubmit = async (e) => {
    e.preventDefault();
    var formElement = document.querySelector('form');
    var formData = new FormData(formElement);
    // Change spaces to underscores for each sequence
    var seq1 = formData.get('seq1').replaceAll(' ', '_');
    var seq2 = formData.get('seq2').replaceAll(' ', '_');
    // Parse score schema into integers
    var match = parseInt(formData.get('match'));
    var mismatch = parseInt(formData.get('mismatch'));
    var gap = parseInt(formData.get('gap'));
    // Run the algorithm selected
    if (document.getElementById('sw').checked) {
        smithWaterman(seq1, seq2, match, mismatch, gap);
    } else {
        needlemanWunsch(seq1, seq2, match, mismatch, gap);
    }
};


// Global variables
class Node {
    constructor() {
        this.score = 0;
        this.pointers = new Array;
    }
}

var allPaths = new Array();


// Smith-Waterman local alignment algorithm
function smithWaterman(seq1, seq2, match, mismatch, gap) {
    seq1 = Array.from(seq1);
    seq2 = Array.from(seq2);
    match = parseInt(match);
    mismatch = parseInt(mismatch);
    gap = parseInt(gap);

    // Initialize a matrix of Nodes
    var matrix = [];
    for (var i = 0; i < seq1.length + 1; i++) {
        matrix.push([])
        for (var j = 0; j < seq2.length + 1; j++) {
            matrix[i].push(new Node)
        }
    }

    // Fill in the rest of the scores
    for (var i = 1; i < seq1.length + 1; i++) {
        for (var j = 1; j < seq2.length + 1; j++) {
            var lookLeft = matrix[i - 1][j].score + gap;
            if (lookLeft < 0) {
                lookLeft = 0;
            }

            var lookUp = matrix[i][j - 1].score + gap;
            if (lookUp < 0) {
                lookUp = 0;
            }

            var lookDiag = matrix[i - 1][j - 1].score;
            if (seq1[i - 1] === seq2[j - 1]) {
                lookDiag += match;
            } else {
                lookDiag += mismatch;
            }
            if (lookDiag < 0) {
                lookDiag = 0;
            }

            var max = Math.max(lookLeft, lookUp, lookDiag);
            matrix[i][j].score = max;

            // Add the pointers
            if (max !== 0) {
                if (max === lookLeft) {
                    matrix[i][j].pointers.push([-1, 0]);
                }
                if (max === lookUp) {
                    matrix[i][j].pointers.push([0, -1]);
                }
                if (max === lookDiag) {
                    matrix[i][j].pointers.push([-1, -1]);
                }
            }
        }
    }

    // Find the max score in the matrix
    var maxScore = 0;
    var coords = new Array();
    for (var i = 1; i < seq1.length + 1; i++) {
        for (var j = 1; j < seq2.length + 1; j++) {
            if (matrix[i][j].score > maxScore) {
                maxScore = matrix[i][j].score;
            }
        }
    }

    // Find the coordinates of each occurrence of the max score in the matrix
    for (var i = 1; i < seq1.length + 1; i++) {
        for (var j = 1; j < seq2.length + 1; j++) {
            if (matrix[i][j].score === maxScore) {
                coords.push([i, j]);
            }
        }
    }

    // Trace-back
    // Find all possible paths from max score coordinate(s)
    allPaths = new Array();
    for (var x = 0; x < coords.length; x++) {
        // Trace back from the highest score until a score of 0
        var start = coords[x];
        var path = [start]
        var i = start[0];
        var j = start[1];
        while (matrix[i][j].score !== 0) {
            var temp = matrix[i][j].pointers[0];
            i = temp[0] + i;
            j = temp[1] + j;
            path.splice(0, 0, [i, j]);
        }
        allPaths.push(path)
    }

    // Dynamically generate HTML output
    buildMatrixEl(seq1, seq2, matrix, allPaths);
    buildAlignmentsEl(seq1, seq2, matrix, allPaths);
};


// Needleman-Wunsch global alignment algorithm
function needlemanWunsch(seq1, seq2, match, mismatch, gap) {
    seq1 = Array.from(seq1);
    seq2 = Array.from(seq2);
    match = parseInt(match);
    mismatch = parseInt(mismatch);
    gap = parseInt(gap);

    // Initialize a matrix of Nodes
    var matrix = [];
    for (var i = 0; i < seq1.length + 1; i++) {
        matrix.push([])
        for (var j = 0; j < seq2.length + 1; j++) {
            matrix[i].push(new Node)
        }
    }

    // Initialize the first row and column with appropriate additive gap score
    for (var i = 1; i < seq1.length + 1; i++) {
        matrix[i][0].score = gap * i;
    }
    for (var j = 1; j < seq2.length + 1; j++) {
        matrix[0][j].score = gap * j;
    }

    // Fill in the rest of the scores
    for (var i = 1; i < seq1.length + 1; i++) {
        for (var j = 1; j < seq2.length + 1; j++) {
            var lookLeft = matrix[i - 1][j].score + gap;

            var lookUp = matrix[i][j - 1].score + gap;

            var lookDiag = matrix[i - 1][j - 1].score;
            if (seq1[i - 1] === seq2[j - 1]) {
                lookDiag += match;
            } else {
                lookDiag += mismatch;
            }

            var max = Math.max(lookLeft, lookUp, lookDiag);
            matrix[i][j].score = max;

            // Add the pointers
            if (max === lookLeft) {
                matrix[i][j].pointers.push([-1, 0]);
            }
            if (max === lookUp) {
                matrix[i][j].pointers.push([0, -1]);
            }
            if (max === lookDiag) {
                matrix[i][j].pointers.push([-1, -1]);
            }

        }
    }

    // Trace-back from the bottom right to the top left of the matrix
    var start = new Array([seq1.length, seq2.length])
    allPaths = new Array()
    traceBack(start, matrix);
    // Reverse each path in allPaths
    for (var x = 0; x < allPaths.length; x++) {
        allPaths[x] = allPaths[x].reverse()
    }

    // Dynamically generate HTML output
    buildMatrixEl(seq1, seq2, matrix, allPaths);
    buildAlignmentsEl(seq1, seq2, matrix, allPaths);
};


// Recursively trace-back through the matrix
function traceBack(paths, matrix) {
    var i = paths[paths.length - 1][0];
    var j = paths[paths.length - 1][1];
    var iNext = i;
    var jNext = j;
    var pointers = matrix[i][j].pointers;

    // Base case when trace-back reaches the top left of matrix
    if (i === 0 && j === 0) {
        return allPaths.push(paths)
    }

    if (pointers.length === 1) {
        iNext += pointers[0][0];
        jNext += pointers[0][1];
        paths[paths.length] = [iNext, jNext];
        traceBack(paths, matrix);
    } else {
        var pathLen = paths.length
        for (var x = 0; x < pointers.length; x++) {
            iNext = i;
            jNext = j;
            iNext += pointers[x][0];
            jNext += pointers[x][1];
            if (x === 0) {
                paths[paths.length] = [iNext, jNext];
                traceBack(paths, matrix);
            } else {
                paths = allPaths[allPaths.length - 1].slice(0, pathLen);
                paths[paths.length] = [iNext, jNext];
                traceBack(paths, matrix);
            }
        }
    }
};


// Generate the alignments element
function buildAlignmentsEl(seq1, seq2, matrix, allPaths) {
    var results = new Array();

    console.log(matrix)

    // Create formatted strings for each alignment from its path in the matrix
    for (var z = 0; z < allPaths.length; z++) {
        var seq1Results = new Array();
        var alignSymbols = new Array();
        var seq2Results = new Array();

        var path = allPaths[z]
        console.log(path)
        for (var i = 1; i < path.length; i++) {
            if (path[i][0] === 0 && path[i][1] === 0) {
                break;
            }
            var x = path[i - 1][0] - path[i][0];
            var y = path[i - 1][1] - path[i][1];
            var pointer = [x, y];
            console.log(pointer)

            if (pointer[0] === -1 && pointer[1] === -1) {
                if (seq1[path[i][0] - 1] == seq2[path[i][1] - 1]) { //Match
                    seq1Results.push(seq1[path[i][0] - 1]);
                    alignSymbols.push("|");
                    seq2Results.push(seq2[path[i][1] - 1]);
                } else {                                            //Mismatch
                    seq1Results.push(seq1[path[i][0] - 1]);
                    alignSymbols.push("*");
                    seq2Results.push(seq2[path[i][1] - 1]);
                };
            } else if (pointer[0] == 0 && pointer[1] == -1) {       //Gap seq1
                seq1Results.push("-");
                alignSymbols.push(" ");
                seq2Results.push([seq2[path[i][1] - 1]]);
            } else {                                                //Gap seq2
                seq1Results.push(seq1[path[i][0] - 1]);
                alignSymbols.push(" ");
                seq2Results.push("-");
            };
        };
        results.push([seq1Results, alignSymbols, seq2Results])
    };

    // Reset the alignments element
    var alignmentsResults = document.getElementById('alignments');
    alignmentsResults.innerHTML = '';
    // Fill in the alignments element
    for (var x = 0; x < allPaths.length; x++) {
        var tr = document.createElement('tr');
        alignmentsResults.appendChild(tr);
        var td = document.createElement('td');
        td.className = 'selectable alignment-'.concat(x);
        td.id = x;
        tr.appendChild(td);
        var code1 = document.createElement('code');
        code1.innerHTML = results[x][0].join('&nbsp;').concat('<br>');
        td.appendChild(code1)
        var code2 = document.createElement('code');
        code2.innerHTML = results[x][1].join('&nbsp;').concat('<br>');
        td.appendChild(code2)
        var code3 = document.createElement('code');
        code3.innerHTML = results[x][2].join('&nbsp;').concat('<br>');
        td.appendChild(code3)
    }

    // Add event listeners to the alignments
    var alignments = document.getElementsByClassName('selectable');
    for (var i = 0; i < alignments.length; i++) {
        alignments[i].addEventListener('click', highlightMatrix);
    }
};


// Highlight the selected alignment in the matrix element
function highlightMatrix() {
    // Get the path that was clicked
    var alignment = this.classList;
    path_n = parseInt(alignment[1].slice(-1));
    var path = allPaths[path_n]

    // Clear any previously selected alignment
    var toClear = document.getElementsByClassName('selectable');
    for (var x = 0; x < toClear.length; x++) {
        if (toClear.item(x).classList.contains('selected')) {
            toClear.item(x).classList.remove('selected');
        }
    }

    // Highlight the selected alignment
    var toSelect = document.getElementsByClassName(alignment[1]);
    toSelect.item(0).classList.add('selected');

    // Clear the matrix of any highlighted cells
    var matrixCells = document.getElementsByClassName('matrix-cell');
    for (var x = 0; x < matrixCells.length; x++) {
        if (matrixCells.item(x).classList.contains('selected')) {
            matrixCells.item(x).classList.remove('selected');
        }
    }

    // Highlight the selected path
    for (var x = 0; x < path.length; x++) {
        var cell = path[x];
        var col = cell[0] + 1;
        var row = cell[1] + 1;
        var matrixItem = document.getElementsByClassName('item-'.concat(col, '-', row));
        matrixItem.item(0).classList.add('selected');
    }

}


// Generate the matrix element
function buildMatrixEl(seq1, seq2, matrix, allPaths) {
    // Reset the matrix element
    matrixContainer = document.getElementById('matrix-container');
    matrixContainer.innerHTML = '';

    // Fill in the matrix element
    for (var i = 0; i < seq1.length + 2; i++) {
        var matrixRow = document.createElement('div');
        matrixRow.className = 'col-'.concat(i);
        matrixContainer.appendChild(matrixRow);

        for (var j = 0; j < seq2.length + 2; j++) {
            var matrixItem = document.createElement('div');
            matrixItem.className = 'matrix-cell item-'.concat(i, '-', j);

            if (i === 0 && j < 2) {                 // First 2 spaces are blank
                matrixItem.className += ' seq-char';
                matrixItem.innerHTML = '&nbsp;';
            } else if (i === 0 && j > 1) {          // Show seq1 across top
                matrixItem.className += ' seq-char';
                matrixItem.innerHTML = seq2[j - 2];
            } else if (i > 0 && j > 0) {            // Fill in matrix scores
                matrixItem.className += ' score';
                matrixItem.innerHTML = matrix[i - 1][j - 1].score;
            } else if (i < 2 && j === 0) {          // First 2 row spaces blank
                matrixItem.className += ' seq-char';
                matrixItem.innerHTML = '&nbsp;';
            } else {                                // Show seq2 down the left
                matrixItem.className += ' seq-char';
                matrixItem.innerHTML = seq1[i - 2];
            }

            matrixRow.appendChild(matrixItem);
        }
    }

    // Variables for dynamically adjusting the matrix in CSS
    var columns = seq1.length + 2;
    var rows = seq2.length + 2;
    document.documentElement.style.setProperty('--colNumMatrix', columns);
    document.documentElement.style.setProperty('--rowNumMatrix', rows);
};
