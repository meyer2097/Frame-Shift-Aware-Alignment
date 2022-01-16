"""
@author: Jules Kreuer, Catalina, Samira Breitling, Sebastian Meyer
"""

from numpy import zeros, int_
import numpy.typing as npt
from os.path import isfile

from .fasta import read_fasta, translate_seq
from .blosum import BLOSUM


def __alignment_core(dnaSeq: str, aaSeq: str,
                     gap: int, shift: int,
                     bm: BLOSUM):
    """
    Core implementation of the frameshift aware Needleman-Wunsch for DNA and AA-Sequences.

    Parameters:
        dnaSeq:  string, dna sequence.
        aaSeq:   string, amino-acid sequence.
        gap:     int, score for gap.
        shift:   int, score for frameshift.
        bm:      BLOSUM Object, blosum matrix for AA-alignment scores.
    Returns:
        score: Int, Score of aligment
        dnaSeq_align: String, alignment of the DNA sequence.
            With \\ denoting a backward-framshift.
            With / denoting a forward-framshift.
            With - denoting a gap.
        aaSeq_align: String, alignment of the AA sequence.
            With - denoting a gap.
    """

    n = len(aaSeq) + 1   # Row
    m = len(dnaSeq) + 3  # Col

    # Add space to compensate index shifting for amino-acid
    aaSeq = f" {aaSeq}"

    # Init matrix with zeros
    score_matrix = zeros((n, m), dtype=int)  # type: npt.NDArray[int_]

    # Init traceback matrix
    traceback_matrix = zeros((n, m), dtype=str)  # type: npt.NDArray[int_]

    # First three columns
    for i in range(n):
        score_matrix[i][0] = -i*gap
        traceback_matrix[i][0] = "U"

        score_matrix[i][1] = -i*gap - shift
        traceback_matrix[i][1] = "1"

        score_matrix[i][2] = -i*gap - shift
        traceback_matrix[i][2] = "2"

    # First row
    for j in range(0, m, 3):
        score_matrix[0][j] = -j*gap
        try:
            score_matrix[0][j+1] = -j*gap - shift
            score_matrix[0][j+2] = -j*gap - shift
        except IndexError:
            pass

        traceback_matrix[0][j] = "3"
        try:
            traceback_matrix[0][j+1] = "1"
            traceback_matrix[0][j+2] = "2"
        except IndexError:
            pass

    # Actual DP
    for i in range(1, n):      # Row, AA
        for j in range(3, m):  # Col, DNA

            translated = translate_seq(dnaSeq[j-3:j])
            align_score = bm.get(translated, aaSeq[i])

            # Regular alignment
            align_dna_aa = (score_matrix[i-1][j-3] + align_score, "D")
            insert_amino = (score_matrix[i-1][j] - gap,           "U")
            l3_dna_insert = (score_matrix[i][j-3] - gap,          "3")

            # Frameshift
            l2_dna_insert = (score_matrix[i][j-2] - shift, "2")
            l1_dna_insert = (score_matrix[i][j-1] - shift, "1")

            # Take best option
            max_path = max(align_dna_aa,
                           insert_amino,
                           l3_dna_insert,
                           l2_dna_insert,
                           l1_dna_insert,
                           key=lambda t: float(t[0]))

            score_matrix[i][j] = int(max_path[0])
            traceback_matrix[i][j] = max_path[1]

    # Traceback
    dnaSeq_align = ""
    aaSeq_align = ""

    i = n-1
    j = m-3

    # Go from bottom right (-3), to top left
    while i > 1 or j > 1:
        # Amino Match / Missmatch
        if traceback_matrix[i][j] == "D":
            dnaSeq_align = translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
            aaSeq_align = aaSeq[i] + aaSeq_align
            j -= 3
            i -= 1

        # Amino insert
        elif traceback_matrix[i][j] == "U":
            dnaSeq_align = "-" + dnaSeq_align
            aaSeq_align = aaSeq[i] + aaSeq_align
            i -= 1

        # L3 insert
        elif traceback_matrix[i][j] == "3":
            dnaSeq_align = translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
            aaSeq_align = "-" + aaSeq_align
            j -= 3

        # L2 insert
        elif traceback_matrix[i][j] == "2":
            dnaSeq_align = "\\" + dnaSeq_align
            aaSeq_align = "-" + aaSeq_align
            j -= 2

        # L1 insert
        elif traceback_matrix[i][j] == "1":
            dnaSeq_align = "/" + dnaSeq_align
            aaSeq_align = "-" + aaSeq_align
            j -= 1

        else:
            Exception("Traceback Traversal Error")

    return score_matrix[n-1][m-3], dnaSeq_align, aaSeq_align


def align(dnaSeq: str, aaSeq: str, gap: int, shift: int,
          blosum, out=False, verbose: bool = False):

    """
    Frameshift aware Needleman-Wunsch for DNA and AA-Sequences.

    Parameters:
        dnaSeq:  string, dna sequence or path to a single entry fasta file.
        aaSeq:   string, amino-acid sequence or path to a single entry fasta file.
        gap:     int, gap penalty (absolute value).
        shift:   int, frameshift penalty (absolute value).
        bm:      BLOSUM Object, blosum matrix for AA-alignment scores.
        verbose: bool, print score and alignment.
    Returns:
        score: Int, Score of aligment
        dnaSeq_align: String, alignment of the DNA sequence.
            With \\ denoting a backward-framshift.
            With / denoting a forward-framshift.
            With - denoting a gap.
        aaSeq_align: String, alignment of the AA sequence.
            With - denoting a gap.
    """

    if isfile(dnaSeq):
        dnaSeq = read_fasta(dnaSeq)[0]

    if isfile(aaSeq):
        aaSeq = read_fasta(aaSeq)[0]

    bm = BLOSUM(blosum)

    if shift <= 0:
        Warning("Shift penalty lower equal zero!")
    if gap <= 0:
        Warning("Gap penalty lower equal zero!")
    if shift < gap:
        Warning("Gap penalty is larger than shift penalty")

    score, dnaSeq_align, aaSeq_align = __alignment_core(dnaSeq, aaSeq, gap, shift, bm)

    if verbose:
        print(f"Score: {score}")
        print(dnaSeq_align)
        print(aaSeq_align)

    # TODO IO
    # PRINT Output
    # WRITE Output

    return score, dnaSeq_align, aaSeq_align
