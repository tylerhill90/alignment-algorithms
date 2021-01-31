seqAlignForm.onsubmit = async (e) => {
    e.preventDefault();
    var formElement = document.querySelector('form');
    var formData = new FormData(formElement);
    var seq1 = formData.get('seq1');
    var seq2 = formData.get('seq2');
    var match = parseInt(formData.get('match'));
    var mismatch = parseInt(formData.get('mismatch'));
    var gap = parseInt(formData.get('gap'));
    if (document.getElementById('sw').checked) {
        smithWaterman(seq1, seq2, match, mismatch, gap);
    } else {
        needlemanWunsch(seq1, seq2, match, mismatch, gap);
    }
};


class Node {
    constructor() {
        this.score = 0;
        this.pointers = new Array;
    }
}


var allPaths = new Array();


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

    // Add the scores
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

    buildMatrixEl(seq1, seq2, matrix, allPaths);

    buildAlignmentsEl(seq1, seq2, matrix, allPaths);
};


function buildAlignmentsEl(seq1, seq2, matrix, allPaths) {
    var results = new Array();

    for (var x = 0; x < allPaths.length; x++) {
        var seq1Results = new Array();
        var alignSymbols = new Array();
        var seq2Results = new Array();

        var path = allPaths[x]

        for (var i = 1; i < path.length; i++) {
            var pointer = matrix[path[i][0]][path[i][1]].pointers[0];
            if (pointer[0] === -1 && pointer[1] === -1) {
                if (seq1[path[i][0] - 1] == seq2[path[i][1] - 1]) { //Match
                    seq1Results.push(seq1[path[i][0] - 1]);
                    alignSymbols.push("|");
                    seq2Results.push(seq2[path[i][1] - 1]);
                } else {                                            //Mismatch
                    seq1Results.push(seq1[path[i][0] - 1]);
                    alignSymbols.push(" ");
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

    var alignments = document.getElementsByClassName('selectable');
    for (var i = 0; i < alignments.length; i++) {
        alignments[i].addEventListener('click', highlightMatrix);
    }
};


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
    for (var x = 1; x < path.length; x++) {
        var cell = path[x];
        var col = cell[0] + 1;
        var row = cell[1] + 1;
        var matrixItem = document.getElementsByClassName('item-'.concat(col, '-', row));
        matrixItem.item(0).classList.add('selected');
    }

}


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

    var columns = seq1.length + 2;
    var rows = seq2.length + 2;
    document.documentElement.style.setProperty('--colNumMatrix', columns);
    document.documentElement.style.setProperty('--rowNumMatrix', rows);
};


function needlemanWunsch(seq1, seq2, match, mismatch, gap) {
    seq1 = Array.from(seq1);
    seq2 = Array.from(seq2);

    // Initialize an empty matrix
    var scoreMatrix = [...Array(seq1.length + 1)].map(e => Array(seq2.length + 1));

    // Fill in the first row and column with the appropriate values
    scoreMatrix[0][0] = 0;
    for (var i = 1; i < seq1.length + 1; i++) {
        scoreMatrix[i][0] = gap * i
    };
    for (var i = 1; i < seq2.length + 1; i++) {
        scoreMatrix[0][i] = gap * i
    };
};
