"""
@author: Jules Kreuer, Catalina, Samira Breitling, Sebastian Meyer
"""

import argparse
from numpy import zeros
from os.path import isfile

import blosum as bl
import miniFasta as ft

def __alignment_core(dnaSeq: str, aaSeq: str,
                     gap_open: int, gap_extend: int, shift: int,
                     bm: bl.BLOSUM):
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

    # TODO: 

    n = len(aaSeq) + 1  # Row
    m = len(dnaSeq) + 3  # Col

    # Add space to compensate index shifting for amino-acid
    aaSeq = f" {aaSeq}"

    # Init matrix with zeros
    score_matrix = zeros((n, m), dtype=int)
    
    top_matrix = zeros((n, m), dtype=float)    #xval
    bottom_matrix = zeros((n, m), dtype=float) #yval

    # Init traceback matrix
    traceback_matrix = zeros((n, m), dtype=str)

    # First three columns
    for i in range(1, n):
        basePenalty = -3*gap_open - (i*gap_extend)
        score_matrix[i][0] = basePenalty
        traceback_matrix[i][0] = "U"
        top_matrix[i][0] = float("-inf")

        score_matrix[i][1] = basePenalty - shift
        traceback_matrix[i][1] = "1"
        top_matrix[i][1] = float("-inf")

        score_matrix[i][2] = basePenalty - shift
        traceback_matrix[i][2] = "2"
        top_matrix[i][2] = float("-inf")


    # First row
    for j in range(0, m, 3):
        basePenalty = -gap_open - (j*gap_extend)

        score_matrix[0][j] = basePenalty

        try:
            score_matrix[0][j+1] = basePenalty - shift
            score_matrix[0][j+2] = basePenalty - shift
        except IndexError:
            pass

        traceback_matrix[0][j] = "3"
        try:
            traceback_matrix[0][j+1] = "1"
            traceback_matrix[0][j+2] = "2"
        except IndexError:
            pass

        bottom_matrix[0][j] = float("-inf")
        try:
            bottom_matrix[0][j+1] = float("-inf")
            bottom_matrix[0][j+2] = float("-inf")
        except IndexError:
            pass

    # Actual DP
    for i in range(1, n):  # Row, AA
        for j in range(3, m):  # Col, DNA

            translated_seq = ft.translate_seq(dnaSeq[j-3:j])
            align_score = bm[f"{translated_seq}{aaSeq[i]}"]

            # Regular alignment
            align_dna_aa = (max(score_matrix[i-1][j-3],
                                top_matrix[i-1][j-3],
                                bottom_matrix[i-1][j-3]) + align_score, "D")
            
            #align_dna_aa = (score_matrix[i-1][j-3] + align_score, "D")

            # Affine Gaps
            insert_amino_open = (score_matrix[i-1][j] - gap_open,  "U")
            insert_amino_extend = (top_matrix[i-1][j] - gap_extend,"U")
            insert_amino = max(insert_amino_open, insert_amino_extend)
        
            l3_dna_insert_open = (score_matrix[i][j-3] - gap_open,      "3")
            l3_dna_insert_extend = (bottom_matrix[i][j-3] - gap_extend, "3")
            l3_dna_insert = max(l3_dna_insert_open, l3_dna_insert_extend)

            # Frameshift
            l2_dna_insert = (score_matrix[i][j-2] - shift, "2")
            l1_dna_insert = (score_matrix[i][j-1] - shift, "1")

            # Take best option
            max_path = max(align_dna_aa,
                           insert_amino,
                           l3_dna_insert,
                           l2_dna_insert,
                           l1_dna_insert)
                           
            score_matrix[i][j]     = int(max_path[0])
            top_matrix[i][j]       = int(insert_amino[0])
            bottom_matrix[i][j]    = int(l3_dna_insert[0])
            traceback_matrix[i][j] = max_path[1]

    # Traceback
    dnaSeq_align = ""
    aaSeq_align = ""

    print(score_matrix)
    print(top_matrix)
    print(bottom_matrix)
    print(traceback_matrix)
    i = n-1
    j = m-3

    # Go from bottom right (-3), to top left
    while i > 1 or j > 1:
        # Amino Match / Missmatch
        if traceback_matrix[i][j] == "D":
            dnaSeq_align = ft.translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
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
            dnaSeq_align = ft.translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
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


def align(dnaSeq: str, aaSeq: str, gap_open: int, gap_extend:int, shift: int,
          blosum_number, out=False, verbose: bool = False):

    """
    Frameshift aware Needleman-Wunsch for DNA and AA-Sequences.

    Parameters:
        dnaSeq:  string, dna sequence.
        aaSeq:   string, amino-acid sequence.
        gap:     int, score for gap.
        shift:   int, score for frameshift.
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
        dnaSeq = ft.read(dnaSeq)[0]

    if isfile(aaSeq):
        aaSeq = ft.read(aaSeq)[0]

    bm = bl.BLOSUM(blosum_number)

    if shift <= 0:
        Warning("Shift penalty lower equal zero!")
    if gap_open <= 0:
        Warning("Gap penalty lower equal zero!")
    if 2*shift < gap_open:
        Warning("Gap penalty is larger than 2*shift penalty")

    score, dnaSeq_align, aaSeq_align = __alignment_core(dnaSeq, aaSeq, gap_open, gap_extend, shift, bm)

    if verbose:
        print(f"Score: {score}")
        print(dnaSeq_align)
        print(aaSeq_align)

    # TODO IO
    # PRINT Output
    # WRITE Output

    return score, dnaSeq_align, aaSeq_align


align("ATGATGATGCCC",
      "PMP",
          gap_open=7,
          gap_extend=1,
          shift=30,
          blosum_number=62,
          out=False,
          verbose=True)

if __name__ == "__main__" and False:

    # Parse arguments
    parser = argparse.ArgumentParser(description='A NICE TEXT')
    parser.add_argument('-d', "--dnaseq", action='store', dest='dnaseq',
                        help="DNA Sequence to path to FASTA file with one entry.",
                        required=True)

    parser.add_argument('-aa', "--aaseq", action='store', dest='aaseq',
                        help="Amino-acid sequence or path to FASTA file with one entry.",
                        required=True)

    parser.add_argument('-gp', "--gap", action='store', dest='gap',
                        help="DNA gap penalty",
                        required=True)

    parser.add_argument('-sp', "--shift", action='store', dest='shift',
                        help="Frameshift penalty",
                        required=True)

    parser.add_argument('-b', "--blosum", action='store', dest='blosum', default=62,
                        help=""" Specify blosum matrix. One of: 45, 50, 62, 80, 90. Default: 62.
                                 Will be irgnored if -bp is set.""",
                        required=False)

    parser.add_argument('-bp', "--blosum_path", action='store', dest='blosum_path',
                        help="Specify path to blosum matrix. Use only if -b is not a viable solution.",
                        required=False)

    parser.add_argument('-o', "--out", action='store', dest='out',
                        help="Specify path to output file.",
                        required=False)

    args = parser.parse_args()

    arg_dnaseq = args.dnaseq
    arg_aaseq = args.aaseq
    arg_gap = int(args.gap)
    arg_shift = int(args.shift)
    arg_out = args.out

    if args.blosum_path:
        arg_blosum = args.blosum_path
    else:
        arg_blosum = int(args.blosum)

    align(arg_dnaseq,
          arg_aaseq,
          arg_gap,
          arg_shift,
          arg_blosum,
          arg_out,
          verbose=True)
