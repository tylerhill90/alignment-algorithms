seqAlignForm.onsubmit = async (e) => {
    e.preventDefault();
    var formElement = document.querySelector('form');
    var formData = new FormData(formElement);
    var seq1 = formData.get('seq1');
    var seq2 = formData.get('seq2');
    var match = formData.get('match');
    var mismatch = formData.get('mismatch');
    var gap = formData.get('gap');
    if (document.getElementById('sw').checked) {
        var alignment = smithWaterman(seq1, seq2, match, mismatch, gap);
    } else {
        var alignment = needlemanWunsch(seq1, seq2, match, mismatch, gap);
    }
};

function smithWaterman(seq1, seq2, match, mismatch, gap) {
    seq1 = Array.from(seq1)
    seq2 = Array.from(seq2)
    match = parseInt(match)
    mismatch = parseInt(mismatch)
    gap = parseInt(gap)

    // Initialize an empty matrix for score and for remembering pointers
    var scoreMatrix = [...Array(seq1.length + 1)].map(() => Array(seq2.length + 1));
    var pointerMatrix = [...Array(seq1.length + 1)].map(() => Array(seq2.length + 1));

    // Fill in the first row and column with zeros
    for (var i = 0; i < seq1.length + 1; i++) {
        scoreMatrix[i][0] = 0
    };
    for (var i = 1; i < seq2.length + 1; i++) {
        scoreMatrix[0][i] = 0
    };

    // Fill in the rest of the score and pointer matrix
    for (var i = 1; i < seq1.length + 1; i++) {
        for (var j = 1; j < seq2.length + 1; j++) {
            var lookLeft = scoreMatrix[i - 1][j] + gap;
            if (lookLeft < 0) {
                lookLeft = 0;
            };

            var lookUp = scoreMatrix[i][j - 1] + gap;
            if (lookUp < 0) {
                lookUp = 0;
            };

            var lookDiag = scoreMatrix[i - 1][j - 1];
            if (seq1[i - 1] === seq2[j - 1]) {
                lookDiag += match;
            } else {
                lookDiag += mismatch;
            };
            if (lookDiag < 0) {
                lookDiag = 0;
            };

            var max = Math.max(lookLeft, lookUp, lookDiag);
            scoreMatrix[i][j] = max;

            // Fill in the pointer matrix
            if (max !== 0) {
                if (max === lookLeft) {
                    pointerMatrix[i][j] = [-1, 0];
                };
                if (max === lookUp) {
                    if (pointerMatrix[i][j]) {
                        pointerMatrix[i][j].push([0, -1])
                    } else {
                        pointerMatrix[i][j] = [0, -1];
                    };
                };
                if (max === lookDiag) {
                    if (pointerMatrix[i][j]) {
                        pointerMatrix[i][j].push([-1, -1])
                    } else {
                        pointerMatrix[i][j] = [-1, -1];
                    };
                };
            };
        };
    };

    // Find the max score and its coordinates in the matrix
    var maxScore = 0;
    var coords = [0, 0]
    for (var i = 1; i < seq1.length + 1; i++) {
        for (var j = 1; j < seq2.length + 1; j++) {
            if (scoreMatrix[i][j] > maxScore) {
                maxScore = scoreMatrix[i][j];
                coords = [i, j]
            }
        }
    }

    // Trace back from the highest score until a score of 0
    var path = [coords]
    var i = coords[0];
    var j = coords[1];
    while (scoreMatrix[coords[0]][coords[1]] !== 0) {
        i = coords[0];
        j = coords[1];
        var iNext = pointerMatrix[i][j][0] + i;
        var jNext = pointerMatrix[i][j][1] + j;
        coords = [iNext, jNext]

        path.splice(0, 0, coords)
    }

    // Print the alignment
    var seq1Results = new Array();
    var alignSymbols = new Array();
    var seq2Results = new Array();

    for (var i = 1; i < path.length; i++) {
        var pointer = pointerMatrix[path[i][0]][path[i][1]];

        if (pointer[0] == -1 && pointer[1] == -1) {
            seq1Results.push(seq1[path[i][0] - 1]);
            alignSymbols.push("|");
            seq2Results.push(seq2[path[i][1] - 1]);
        } else if (pointer[0] == 0 && pointer[1] == -1) {
            seq1Results.push("-");
            alignSymbols.push(" ");
            seq2Results.push([seq2[path[i][1] - 1]]);
        } else {
            seq1Results.push(seq1[path[i][0] - 1]);
            alignSymbols.push("|");
            seq2Results.push("-");
        };
    };

    var alignmentEl = document.getElementById('alignment');
    alignmentEl.innerHTML = ""
    var n = seq1Results.length;
    for (var i = 0; i < n; i++) {
        let gridItem = document.createElement('div');
        gridItem.className = 'grid-item';
        gridItem.innerHTML += seq1Results[i];
        alignmentEl.appendChild(gridItem);
    };

    for (var i = 0; i < n; i++) {
        let gridItem = document.createElement('div');
        gridItem.className = 'grid-item';
        gridItem.innerHTML += alignSymbols[i];
        alignmentEl.appendChild(gridItem);
    };

    for (var i = 0; i < n; i++) {
        let gridItem = document.createElement('div');
        gridItem.className = 'grid-item';
        gridItem.innerHTML += seq2Results[i];
        alignmentEl.appendChild(gridItem);
    };

    document.documentElement.style.setProperty("--colNum", n);

};

function needlemanWunsch(seq1, seq2, match, mismatch, gap) {
    seq1 = Array.from(seq1)
    seq2 = Array.from(seq2)

    // Initialize an empty matrix
    var scoreMatrix = [...Array(seq1.length + 1)].map(e => Array(seq2.length + 1));

    // Fill in the first row and column with the appropriate values
    scoreMatrix[0][0] = 0
    for (var i = 1; i < seq1.length + 1; i++) {
        scoreMatrix[i][0] = gap * i
    };
    for (var i = 1; i < seq2.length + 1; i++) {
        scoreMatrix[0][i] = gap * i
    };
    console.log(scoreMatrix)
};
