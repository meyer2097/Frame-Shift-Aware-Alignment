"""
DP Implementation of a frameshift aware alignment
Version: Catalina & Jules
"""

import numpy as np
import sys
import fasta as ft
import blosum as bl

SCORES = {'gap': 1, 'frame1': 2, 'frame2': 2, 'frame3': 1}


def fa_nw(dnaSeq, aaSeq, bm, verbose=True):
    """
    Frameshift aware implementation of Needleman-Wunsch for DNA and AA-Sequences
    Parameters:
        todo
    Returns:
        todo
    """

    n = len(aaSeq) + 1 # Row
    m = len(dnaSeq)+ 3 # Col

    aaSeq  = f"*{aaSeq}"
    dbaSeq = f"*{dnaSeq}"

    # Init matrix with zeros
    score_matrix     = np.zeros((n, m), dtype = int)

    # Init traceback matrix
    traceback_matrix = np.zeros((n, m), dtype = str)
    
    # First three columns
    for i in range(n):
        score_matrix[i][0] = -i*SCORES['gap']
        traceback_matrix[i][0] = "U"

        score_matrix[i][1] = -i*SCORES['gap'] - SCORES['frame1']
        traceback_matrix[i][1] = "1"

        score_matrix[i][2] = -i*SCORES['gap'] - SCORES['frame2']
        traceback_matrix[i][2] = "2"
        

    # First row 
    for j in range(0, m, 3):
        score_matrix[0][j] = -j*SCORES['frame3']
        try:
            score_matrix[0][j+1] = -j*SCORES['frame3'] - SCORES['frame1']
            score_matrix[0][j+2] = -j*SCORES['frame3'] - SCORES['frame2']
        except IndexError:
            pass

        traceback_matrix[0][j] = "3"
        try:
            traceback_matrix[0][j+1] = "1"
            traceback_matrix[0][j+2] = "2"
        except IndexError:
            pass

    print(score_matrix)
    print(traceback_matrix)

    # Actual DP
    for i in range(1, n): # Row, AA
        for j in range(3, m): # Col, DNA
                
            translated_seq = ft.translate_seq(dnaSeq[j-3:j])
            align_score = bm.get(translated_seq, aaSeq[i])

            # Regular alignment
            align_dna_aa  = (score_matrix[i-1][j-3] + align_score     , "D")
            insert_amino  = (score_matrix[i-1][j]   - SCORES['gap']   , "U")
            l3_dna_insert = (score_matrix[i][j-3]   - SCORES['frame3'], "3")

            # Frameshift
            l2_dna_insert = (score_matrix[i][j-2]   - SCORES['frame2'], "2")
            l1_dna_insert = (score_matrix[i][j-1]   - SCORES['frame1'], "1")
            
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
    
    print(i,j)
    # Go from bottom right, to top left
    while i > 1 or j > 1:
        # Amino Match / Missmatch
        if traceback_matrix[i][j]  == "D":
            dnaSeq_align = ft.translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
            aaSeq_align  = aaSeq[i] + aaSeq_align
            j -= 3
            i -= 1

        # Amino insert
        elif traceback_matrix[i][j] == "U":
            dnaSeq_align = "-" + dnaSeq_align
            aaSeq_align  = aaSeq[i]       + aaSeq_align
            i -= 1
        
        # L3 insert
        elif traceback_matrix[i][j] == "3":
            dnaSeq_align = ft.translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
            aaSeq_align  = "-"  + aaSeq_align
            j -= 3
        
        # L2 insert
        elif traceback_matrix[i][j] == "2":
            dnaSeq_align = "#" + dnaSeq_align
            aaSeq_align  = "-" + aaSeq_align
            j -= 2
        
        # L1 insert
        elif traceback_matrix[i][j] == "1":
            dnaSeq_align = "$"  + dnaSeq_align
            aaSeq_align  = "-"       + aaSeq_align
            j -= 1

        else:
            Exception("Traceback Traversal Error")
    
    if verbose:
        print("Score: " + str(score_matrix[n-1][m-3]))
        print(dnaSeq_align)
        print(aaSeq_align)

    return [dnaSeq_align, aaSeq_align, score_matrix[n-1][m-3]]

bm = bl.BLOSUM(62)
fa_nw("GCCCGCCGCCGCC", "AAAA", bm)
