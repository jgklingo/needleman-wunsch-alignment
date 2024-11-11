class SparseMatrix:
    # i is the row (seq1), j is the column (seq2)
    def __init__(self, seq1length, seq2length):
        self.data = {}  # {(i, j): value, previous_i, previous_j}
        self.ROWS = seq1length + 1
        self.COLS = seq2length + 1
    def getValue(self, index):
        item = self.data.get(index, float('inf'))
        return item[0] if item != float('inf') else float('inf')
    def getPrevious(self, index):
        item = self.data.get(index, None)
        return (item[1], item[2]) if item != None else None
    def setValuePrevious(self, index, value, previous):
        assert len(previous) == 2 and type(previous) == tuple
        i, j = index[0], index[1]
        if not 0 <= i < self.ROWS or not 0 <= j < self.COLS:
            return  #raise IndexError("Index out of bounds")
        self.data[index] = (value, previous[0], previous[1])


def align(
        seq1: str,
        seq2: str,
        match_award=-3,
        indel_penalty=5,
        sub_penalty=1,
        banded_width=-1,
        gap='-'
    ) -> tuple[float, str | None, str | None]:
    """
        Align seq1 against seq2 using Needleman-Wunsch
        Put seq1 on left (j) and seq2 on top (i)
        => matrix[i][j]
        :param seq1: the first sequence to align; should be on the "left" of the matrix
        :param seq2: the second sequence to align; should be on the "top" of the matrix
        :param match_award: how many points to award a match
        :param indel_penalty: how many points to award a gap in either sequence
        :param sub_penalty: how many points to award a substitution
        :param banded_width: banded_width * 2 + 1 is the width of the banded alignment; -1 indicates full alignment
        :param gap: the character to use to represent gaps in the alignment strings
        :return: alignment cost, alignment 1, alignment 2
    """
    def score(data, index):
        if index == (0, 0):
            return 0, (None, None)

        i, j = index[0], index[1]
        if i < 1 or j < 1:
            diff = float('inf')
        elif seq1[i - 1] == seq2[j - 1]:
            diff = match_award
        else:
            diff = sub_penalty
        directions = {
            (-1, -1): data.getValue((i - 1, j - 1)) + diff,  # diagonal
            (0, -1): data.getValue((i, j - 1)) + indel_penalty,  # left
            (-1, 0): data.getValue((i - 1, j)) + indel_penalty  # top
        }
        minimum = float('inf')
        previous = (None, None)
        for key in directions.keys():
            if directions[key] < minimum:
                minimum = directions[key]
                previous = (i + key[0], j + key[1])
        assert previous != (None, None)
        return minimum, previous

    if banded_width == -1:
        banded_width = max(len(seq1) + 1, len(seq2) + 1)
    data = SparseMatrix(len(seq1), len(seq2))
    for diagonal in range(min(len(seq1) + 1, len(seq2) + 1)):
        for k in range(banded_width + 1):
            if 0 <= diagonal + k < data.ROWS and 0 <= diagonal < data.COLS:
                down_score = score(data, (diagonal + k, diagonal))
                data.setValuePrevious((diagonal + k, diagonal), down_score[0], down_score[1])
            if 0 <= diagonal < data.ROWS and 0 <= diagonal + k < data.COLS:
                right_score = score(data, (diagonal, diagonal + k))
                data.setValuePrevious((diagonal, diagonal + k), right_score[0], right_score[1])

    # recover the character-by-character alignment
    alignment1, alignment2 = [], []
    curr = (len(seq1), len(seq2))
    while data.getPrevious(curr) != (None, None):
        prev = curr
        curr = data.getPrevious(curr)
        if curr[0] < prev[0] and curr[1] < prev[1]:
            alignment1.append(seq1[prev[0] - 1])
            alignment2.append(seq2[prev[1] - 1])
        elif curr[0] < prev[0]:
            alignment1.append(seq1[prev[0] - 1])
            alignment2.append(gap)
        elif curr[1] < prev[1]:
            alignment1.append(gap)
            alignment2.append(seq2[prev[1] - 1])
    alignment1, alignment2 = alignment1[::-1], alignment2[::-1]
    
    return (data.getValue((len(seq1), len(seq2))), ''.join(alignment1), ''.join(alignment2))

