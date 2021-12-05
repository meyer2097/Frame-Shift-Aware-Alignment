"""
DP Implementation of a frameshift aware alignment
Version: Catalina & Jules
"""

import numpy as np
import sys
import fasta as ft
import blosum as bl

SCORES = {'gap': 4, 'frame1': 12, 'frame2': 13, 'frame3': 14}


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

        ''' 
        score_matrix[i][1] = -i*SCORES['gap'] - SCORES['frame2']
        traceback_matrix[i][1] = "2"

        score_matrix[i][2] = -i*SCORES['gap'] -SCORES['frame1']
        traceback_matrix[i][2] = "1"
        '''

    # First row 
    for j in range(m):
        score_matrix[0][j] = -j*SCORES['frame3']
        traceback_matrix[0][j] = "3"
    

    print(score_matrix)
    print(traceback_matrix)

    # Actual DP
    for i in range(1, n): # Row, AA
        for j in range(3, m): # Col, DNA
                
            translated_seq = ft.translate_seq(dnaSeq[j-3:j])
            align_score = bm.get(translated_seq, aaSeq[i])

            # Regular alignment
            align_dna_aa  = score_matrix[i-1][j-3] + align_score
            insert_amino  = score_matrix[i-1][j]   - SCORES['gap']  
            l3_dna_insert = score_matrix[i][j-3]   - SCORES['frame3'] 

            # Frameshift
            l2_dna_insert = score_matrix[i][j-2]   - SCORES['frame2']
            l1_dna_insert = score_matrix[i][j-1]   - SCORES['frame1']
            
            # Take best option
            score_matrix[i][j] = max(align_dna_aa, 
                                     insert_amino,
                                     l3_dna_insert,
                                     l2_dna_insert,
                                     l1_dna_insert)

            # Insert path to traceback matrix
            if score_matrix[i][j] == align_dna_aa:
                traceback_matrix[i][j] = "D"
            elif score_matrix[i][j] == insert_amino:
                traceback_matrix[i][j] = "U"
            elif score_matrix[i][j] == l3_dna_insert:
                traceback_matrix[i][j] = "3"
            elif score_matrix[i][j] == l2_dna_insert:
                traceback_matrix[i][j] = "2"
            elif score_matrix[i][j] == l1_dna_insert:
                traceback_matrix[i][j] = "1"
            else:
                # Can never happen
                exit("Traceback Insert Error")
    if verbose:
        print(score_matrix)
        print(traceback_matrix)

    # Traceback
    dnaSeq_align = ""
    aaSeq_align = ""
    i = n-1
    j = m-3
    # Go from bottom right, to top left
    while i > 1 or j > 1:
        # Amino Match / Missmatch
        
        if traceback_matrix[i][j] == "D":
            dnaSeq_align = ft.translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
            aaSeq_align  = aaSeq[i] + aaSeq_align
            j -= 3
            i -= 1

        # Amino Insert
        elif traceback_matrix[i][j] == "U":
            dnaSeq_align = "-" + dnaSeq_align
            aaSeq_align  = aaSeq[i]       + aaSeq_align
            i -= 1
        
        # L3 insert
        elif traceback_matrix[i][j] == "3":
            dnaSeq_align = ft.translate_seq(dnaSeq[j-3:j]) + dnaSeq_align
            aaSeq_align  = "-"  + aaSeq_align
            j -= 3
        
        elif traceback_matrix[i][j] == "2":
            dnaSeq_align = "#" + dnaSeq_align
            aaSeq_align  = "-"       + aaSeq_align
            j -= 2
        
        elif traceback_matrix[i][j] == "1":
            dnaSeq_align = "$"  + dnaSeq_align
            aaSeq_align  = "-"       + aaSeq_align
            j -= 1
        else:
            print(i,j)
            exit("Traceback Traversal Error")
        
        print(i,j)
        print(dnaSeq_align)
        print(aaSeq_align)
    if verbose:
        print("Score: " + str(score_matrix[n-1][m-3]))
        print(dnaSeq_align)
        print(aaSeq_align)

    #return [seq1_top, seq2_bottom, score_matrix[n-1][m-1]]

#nw_regular("AGTAGTAGT", "ACTACTACCCTTTA")
bm = bl.BLOSUM(62)
fa_nw("GCCGCCGCCGCC", "AAAA", bm)
