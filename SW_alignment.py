#!/usr/bin/env python3

# ------------------------------------------------------------------------------
#
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
            local_max = max(matrix[row][column])
            matrix[row][column] = [local_max, local_max, local_max]
            if matrix[row][column][0] < 0:
                matrix[row][column] = [0, 0, 0]

    return matrix


def traverse(o, tree_types=(list, tuple)):
    """Unpack every element in a convuluted list of lists.
    """

    if isinstance(o, tree_types):
        for value in o:
            for subvalue in traverse(value, tree_types):
                yield subvalue
    else:
        yield o


def trace_back(matrix):
    """
    """

    max_cell = max(traverse(matrix))
    start = None
    max_y = -1
    for row in matrix:
        max_x = -1
        max_y += 1
        for cell in row:
            max_x += 1
            if cell[0] == max_cell:
                start = [[max_y, max_x]]

    path = find_path(matrix, start)

    return path


def find_path(matrix, path_taken):
    """
    """

    y, x = path_taken[-1][0], path_taken[-1][1]
    cell_value = matrix[y][x][0]

    # Base case
    if cell_value == 0:
        return path_taken[::-1]

    # Recurse
    up = matrix[y - 1][x][0]
    diag = matrix[y - 1][x - 1][0]
    left = matrix[y][x - 1][0]

    if 0 in [up, diag, left]:
        path_taken.append([0, 0])
        return find_path(matrix, path_taken)
    elif up == max(up, diag, left):
        path_taken.append([y - 1, x])
        return find_path(matrix, path_taken)
    elif diag == max(up, diag, left):
        path_taken.append([y - 1, x - 1])
        return find_path(matrix, path_taken)
    else:
        path_taken.append([y, x - 1])
        return find_path(matrix, path_taken)


def calc_score(seq_1, seq_2, scoring, path):
    """Calculate scores for each alignment given a scoring schema.
    """

    # Assign scoring values
    gap_penalty = scoring[0]
    mismatch_score = scoring[1]
    match_score = scoring[2]

    # Initiate empty strings to store alignment strings
    align_1, match_mismatch, align_2 = "", "", ""

    # Initiate score counter
    score = 0

    # Iterate through each cell except the last
    for cell in range(len(path) - 1):
        # Extract x and y coordinates
        y, x = path[cell][0], path[cell][1]
        # Extract the next x and y coordinates
        y_next, x_next = path[cell + 1][0], path[cell + 1][1]

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

    return str_aligned, score


if __name__ == "__main__":

    # Get scoring scheme from user if not using suggested schema
    while True:
        use_suggested = input("Use suggested scoring schema? y/n\n")
        if use_suggested not in ["Y", "y", "N", "n"]:
            print("Unsupported answer. Try again.")
        else:
            break
    if use_suggested in ["Y", "y"]:
        print("Gap penalty: -2\nMismatch Penalty: -3\nMatch Score: 3")
        scoring = (-2, -3, 3)
    else:
        scoring = define_scoring()

    # Get sequences from user
    sequences = get_sequences()
    seq_1, seq_2 = sequences[0], sequences[1]

    # Initiate the matrix
    matrix = init_matrix(seq_1, seq_2, scoring)

    path = trace_back(matrix)

    alignment = calc_score(seq_1, seq_2, scoring, path)

    print("\nAlignment:")
    print(alignment[0])
    print(f"\nScore: {alignment[1]}")
