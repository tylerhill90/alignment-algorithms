#!/usr/bin/env python3

# ------------------------------------------------------------------------------
# My attempt at the Needleman–Wunsch Algorithm for Global Alignment.
#
# User must input scores for a gap penalty, mismatch, and match as well as two
# sequences to be aligned.
# According to a user's scoring schema this program returns the best
# alignment(s).
#
# All the work is original except for the function traverse() which was taken
# from https://gist.github.com/aljiwala/c77a01f382f5bbc10d2d2b97a7ed0f0a
#
# My work was based on chapter 9 of the book
# "Basic Applied Bioinformatics, First Edition. Chandra Sekhar Mukhopadhyay"
# ------------------------------------------------------------------------------


def define_scoring():
    """Get user input on the scoring schema to be used.
    Input to recieve:
        Gap penalty
        Mismatch score
        Match score
    """

    while True:
        try:
            gap_penalty = input("Gap penalty: ")
            test_gap = int(gap_penalty)
            break
        except ValueError:
            print("Gap penalty must be an integer.")

    while True:
        try:
            mismatch = input("Mismatch penalty: ")
            test_mismatch = int(mismatch)
            break
        except ValueError:
            print("Mismatch penalty must be an integer.")

    while True:
        try:
            match_score = input("Match score: ")
            test_match = int(match_score)
            break
        except ValueError:
            print("Match score must be an integer.")

    return (int(gap_penalty), int(mismatch), int(match_score))


def get_sequences():
    """Get user input for the two sequences to be aligned.
    """

    while True:
        seq_1 = input("1st sequence: ")
        if not seq_1.strip():
            print("Sequence must not be emtpy.")
        else:
            break

    while True:
        seq_2 = input("2nd sequence: ")
        if not seq_2.strip():
            print("Sequence must not be emtpy.")
        else:
            break

    return [seq_1, seq_2]


def init_matrix(seq_1, seq_2, scoring):
    """Construct the alignment matrix based on the sequences provided and
    scoring schema.
    """

    # Define scoring
    gap_penalty, mismatch, match_score = scoring

    # Init a blank matrix with three values for each cell set to 0
    matrix = [[[0, 0, 0] for x in range(len(seq_1) + 1)]
              for x in range(len(seq_2) + 1)]

    # Update matrix's first row and column to include the gap penalty
    # First row
    counter = 0
    for cell in matrix[0]:
        for x in range(len(cell)):
            cell[x] = gap_penalty * counter
        counter += 1

    # First column
    counter = 0
    for line in matrix:
        for first_cell in line:
            for x in range(len(first_cell)):
                first_cell[x] = gap_penalty * counter
            break
        counter += 1

    # Fill in the rest of the matrix
    for row in range(1, len(seq_2) + 1):
        for column in range(1, len(seq_1) + 1):
            for cell in range(3):
                if cell == 0:       # Gap penalty
                    matrix[row][column][cell] = max(
                        matrix[row - 1][column]) + gap_penalty
                elif cell == 1:     # Check if match of mismatch
                    if seq_2[row - 1] == seq_1[column - 1]:
                        matrix[row][column][cell] = max(
                            matrix[row - 1][column - 1]) + match_score
                    else:
                        matrix[row][column][cell] = max(
                            matrix[row - 1][column - 1]) + mismatch
                else:               # Gap penalty
                    matrix[row][column][cell] = max(
                        matrix[row][column - 1]) + gap_penalty

    return matrix


def trace_back(matrix):
    """Find every path through the matrix such that each path starts in the
    bottom right of the matrix and then recursively follows the path towards
    the next cell such that the next cell to follow is indicated by the max
    value(s) for the current cell.
    """

    columns = len(matrix[0]) - 1
    rows = -1
    for row in matrix:
        rows += 1
    start = [[rows, columns]]

    paths = traverse(find_paths(matrix, start))

    paths = [x for x in paths]
    paths = [paths[i:i + 2] for i in range(0, len(paths), 2)]

    new_path = [0]
    for x in range(len(paths)):
        if paths[x] == [0, 0] and x - 1 > 0:
            new_path.append(x)
    new_path.append(len(paths))

    paths = [paths[new_path[i]:new_path[i + 1]]
             for i in range(0, len(new_path) - 1)]

    return paths


def find_paths(matrix, path_taken):
    """Recursively take paths through the matrix as indicated in the trace_back
    function.
    """

    y, x = path_taken[-1][0], path_taken[-1][1]
    current_cell = matrix[y][x]

    if path_taken[-1] == [0, 0]:    # Base case
        return path_taken[::-1]

    else:
        max_value = max(matrix[y][x])

        # Move left or up due to hitting border
        if x == 0:
            return find_paths(matrix, path_taken + [[y - 1, x]])
        elif y == 0:
            return find_paths(matrix, path_taken + [[y, x - 1]])

        # Move all three directions
        elif (max_value == current_cell[0] and max_value == current_cell[1] and
              max_value == current_cell[2]):
            if y - 1 >= 0 and x - 1 >= 0:
                return [
                    find_paths(matrix, path_taken + [[y - 1, x]]),
                    find_paths(matrix, path_taken + [[y, x - 1]]),
                    find_paths(matrix, path_taken + [[y - 1, x - 1]])
                ]

        # Move two directions
        elif max_value == current_cell[0] and max_value == current_cell[1]:
            if y - 1 >= 0 and x - 1 >= 0:
                return [
                    find_paths(matrix, path_taken + [[y - 1, x]]),
                    find_paths(matrix, path_taken + [[y - 1, x - 1]])
                ]
        elif max_value == current_cell[1] and max_value == current_cell[2]:
            if y - 1 >= 0 and x - 1 >= 0:
                return [
                    find_paths(matrix, path_taken + [[y - 1, x - 1]]),
                    find_paths(matrix, path_taken + [[y, x - 1]])
                ]
        elif max_value == current_cell[0] and max_value == current_cell[2]:
            if y - 1 >= 0 and x - 1 >= 0:
                return [
                    find_paths(matrix, path_taken + [[y - 1, x]]),
                    find_paths(matrix, path_taken + [[y, x - 1]])
                ]

        # Move one direction
        elif max_value == current_cell[0]:    # Go up
            if y - 1 >= 0:
                return find_paths(matrix, path_taken + [[y - 1, x]])
        elif max_value == current_cell[1]:  # Go diagonal
            if y - 1 >= 0 and x - 1 >= 0:
                return find_paths(matrix, path_taken + [[y - 1, x - 1]])
        elif max_value == current_cell[2]:  # Go left
            if x - 1 >= 0:
                return find_paths(matrix, path_taken + [[y, x - 1]])


def traverse(o, tree_types=(list, tuple)):
    """Unpack every element in a convuluted list of lists.
    """

    if isinstance(o, tree_types):
        for value in o:
            for subvalue in traverse(value, tree_types):
                yield subvalue
    else:
        yield o


def calc_scores(seq_1, seq_2, scoring, alignments):
    """Calculate scores for each alignment given a scoring schema.
    """

    # Assign scoring values
    gap_penalty = scoring[0]
    mismatch_score = scoring[1]
    match_score = scoring[2]

    # Initiate a scores dictionary to store alignment string and score integer
    scores = {}

    # Iterate through each alignment
    for i in range(len(alignments)):
        # Initiate empty strings to store alignment strings
        align_1, match_mismatch, align_2 = "", "", ""

        # Initiate score counter
        score = 0

        # Iterate through each cell except the last
        for cell in range(len(alignments[i]) - 1):
            # Extract x and y coordinates
            y, x = alignments[i][cell][0], alignments[i][cell][1]
            # Extract the next x and y coordinates
            y_next, x_next = alignments[i][cell +
                                           1][0], alignments[i][cell + 1][1]

            # Determine what direction the path is taking.
            # Update alignment strings and score appropriately.
            if y == y_next and x != x_next:     # Moving right, means gap
                align_1 += seq_1[x_next - 1]
                match_mismatch += " "
                align_2 += "-"
                score += gap_penalty
            elif y != y_next and x == x_next:   # Moving down, means gap
                align_1 += "-"
                match_mismatch += " "
                align_2 += seq_2[y_next - 1]
                score += gap_penalty
            else:                    # Moving diagonal, means match or mismatch
                # Check if sequences match
                if seq_1[x_next - 1] == seq_2[y_next - 1]:
                    align_1 += seq_1[x_next - 1]
                    match_mismatch += "|"
                    align_2 += seq_2[y_next - 1]
                    score += match_score
                # Else they are a mismatch
                else:
                    align_1 += seq_1[x_next - 1]
                    match_mismatch += " "
                    align_2 += seq_2[y_next - 1]
                    score += mismatch_score

        # Format the alignment string
        str_aligned = align_1 + "\n" + match_mismatch + "\n" + align_2

        # Update the scores dictionary
        scores[str_aligned] = score

    return scores


def print_scores(scores):
    """Print the alignment(s) and scores with the max value.
    """

    max_score = max(scores.values())
    counter = 0
    for key, value in scores.items():
        if value == max_score:
            counter += 1
            print(f"\nAlignment {counter}\n" + key)
    print("\nAlignment(s) score: " + str(max_score))


def bias_scores(scores):
    """Bias the scores towards alignments with longer continuous sequences.
    """

    import re

    seq_len = int((len(list(scores.keys())[0]) - 2) / 3)

    for key, value in scores.items():
        alignment = key[seq_len + 1:seq_len * 2 + 1]
        alignment = re.findall(r"\|+", alignment)
        max_chunk = 0
        for chunk in alignment:
            if len(chunk) > max_chunk:
                max_chunk = len(chunk)
        n_chunks = len(alignment)
        if max_chunk != 0:
            scores[key] = value + max_chunk / n_chunks

    return scores


if __name__ == "__main__":

    # Get scoring scheme from user if not using suggested schema
    while True:
        use_suggested = input("Use suggested scoring schema? y/n\n")
        if use_suggested not in ["Y", "y", "N", "n"]:
            print("Unsupported answer. Try again.")
        else:
            break
    if use_suggested in ["Y", "y"]:
        print("Gap penalty: -2\nMismatch Penalty: -1\nMatch Score: 1")
        scoring = (-2, -1, 1)
    else:
        scoring = define_scoring()

    # Get sequences from user
    sequences = get_sequences()
    seq_1, seq_2 = sequences[0], sequences[1]

    # Initiate the matrix
    matrix = init_matrix(seq_1, seq_2, scoring)

    # Print the matrix if it isn't too wide for the screen
    # Uncomment to do so
    if len(matrix[0]) < 10:
        print("\nScoring matrix:")
        for line in matrix:
            for cell in line:
                if len(str(cell)) < 15:
                    print(cell, "\t\t", end="")
                else:
                    print(cell, "\t", end="")
            print()
        print()
    else:
        print("\nScoring matrix too wide to print.\n")

    # Generate alignments
    alignments = trace_back(matrix)

    # Generate scores
    scores = calc_scores(seq_1, seq_2, scoring, alignments)

    # Print normal scores
    print_scores(scores)

    # Bias scores
    biased_scores = bias_scores(scores)

    print("\n\nBIASED ALIGNMENT(S):")
    print("Biased to the longest continuous matching alignment(s)")
    print_scores(biased_scores)
