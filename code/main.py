"""
@author: Jules Kreuer, Catalina, Samira Breitling, Sebastian Meyer
"""

import argparse
import blosum as bl
import fasta  as ft
import numpy as np
from os.path import isfile
def __frameshift_aware_alignment_core(dnaSeq: str, aaSeq: str,
                                    gap: int, shift: int,
                                    bm: bl.BLOSUM(), verbose=False):
    """
    Frameshift aware implementation of Needleman-Wunsch for DNA and AA-Sequences.
    This is the core funcion. 

    Parameters:
        dnaSeq:  string, dna sequence.
        aaSeq:   string, amino-acid sequence.
        gap:     int, score for gap.
        shift:   int, score for frameshift.
        bm:      BLOSUM Object, blosum matrix for AA-alignment scores.
        verbose: bool, default=False, print additional information to stdout.
    Returns:
        score: Int, Score of aligment
        dnaSeq_align: String, alignment of the DNA sequence.
            With \ denoting a backward-framshift.
            With / denoting a forward-framshift.
            With - denoting a gap. 
        aaSeq_align: String, alignment of the AA sequence.
            With - denoting a gap. 
    """

    n = len(aaSeq) + 1 # Row
    m = len(dnaSeq)+ 3 # Col

    # Add space to compensate index shifting
    aaSeq  = f" {aaSeq}"
    dbaSeq = f" {dnaSeq}"

    # Init matrix with zeros
    score_matrix     = np.zeros((n, m), dtype = int)

    # Init traceback matrix
    traceback_matrix = np.zeros((n, m), dtype = str)
    
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
    for i in range(1, n): # Row, AA
        for j in range(3, m): # Col, DNA
                
            translated_seq = ft.translate_seq(dnaSeq[j-3:j])
            align_score = bm.get(translated_seq, aaSeq[i])

            # Regular alignment
            align_dna_aa  = (score_matrix[i-1][j-3] + align_score, "D")
            insert_amino  = (score_matrix[i-1][j]   - gap,         "U")
            l3_dna_insert = (score_matrix[i][j-3]   - gap,         "3")

            # Frameshift
            l2_dna_insert = (score_matrix[i][j-2]   - shift, "2")
            l1_dna_insert = (score_matrix[i][j-1]   - shift, "1")
            
            # Take best option
            max_path = max(align_dna_aa, 
                           insert_amino,
                           l3_dna_insert,
                           l2_dna_insert,
                           l1_dna_insert,
                           key = lambda t: t[0])

            score_matrix[i][j]     = max_path[0]
            traceback_matrix[i][j] = max_path[1]

    if verbose:
        print(score_matrix)
        print(traceback_matrix)

    # Traceback
    dnaSeq_align = ""
    aaSeq_align = ""

    i = n-1
    j = m-3
    
    # Go from bottom right (-3), to top left
    while i > 1 or j > 1:
        # Amino Match / Missmatch
        if traceback_matrix[i][j]  == "D":
            dnaSeq_align = ft.translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
            aaSeq_align  = aaSeq[i] + aaSeq_align
            j -= 3
            i -= 1

        # Amino insert
        elif traceback_matrix[i][j] == "U":
            dnaSeq_align = "-"      + dnaSeq_align
            aaSeq_align  = aaSeq[i] + aaSeq_align
            i -= 1
        
        # L3 insert
        elif traceback_matrix[i][j] == "3":
            dnaSeq_align = ft.translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
            aaSeq_align  = "-" + aaSeq_align
            j -= 3
        
        # L2 insert
        elif traceback_matrix[i][j] == "2":
            dnaSeq_align = "\\" + dnaSeq_align
            aaSeq_align  = "-"  + aaSeq_align
            j -= 2
        
        # L1 insert
        elif traceback_matrix[i][j] == "1":
            dnaSeq_align = "/" + dnaSeq_align
            aaSeq_align  = "-" + aaSeq_align
            j -= 1

        else:
            Exception("Traceback Traversal Error")
    
    return score_matrix[n-1][m-3], dnaSeq_align, aaSeq_align

def frameshift_aware_alignment(dnaSeq: str, aaSeq: str,
                               gap: int, shift: int,
                               blosum, out, verbose: bool):

    
    if isfile(dnaSeq):
        dnaSeq = ft.read_fasta(dnaSeq)[0]
    
    if isfile(aaSeq):
        aaSeq = ft.read_fasta(aaSeq)[0]

    bm = bl.BLOSUM(blosum)
    
    if shift <= 0:
        Warning("Shift penalty lower equal zero!")
    if gap <= 0:
        Warning("Gap penalty lower equal zero!")
    if shift < gap:
        Warning("Gap penalty is larger than shift penalty")

    score, dnaSeq_align, aaSeq_align = __frameshift_aware_alignment_core(dnaSeq, aaSeq, gap, shift, bm, verbose)

    print(f"Score: {score}")
    print(dnaSeq_align)
    print(aaSeq_align)

    # TODO IO
        # PRINT Output
        # WRITE Output

    return

if __name__ == "__main__":

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

    parser.add_argument('-v', "--verbose", action='store_true', dest='verbose', 
        help='Will supress any output.',
        required=False)

    args = parser.parse_args()

    arg_dnaseq  = args.dnaseq
    arg_aaseq   = args.aaseq
    arg_gap     = int(args.gap)
    arg_shift   = int(args.shift)
    arg_out     = args.out
    arg_verbose = args.verbose

    if args.blosum_path:
        arg_blosum = args.blosum_path
    else:
        arg_blosum = int(args.blosum)


    frameshift_aware_alignment(arg_dnaseq,
                               arg_aaseq,
                               arg_gap,
                               arg_shift,
                               arg_blosum,
                               arg_out,
                               arg_verbose)