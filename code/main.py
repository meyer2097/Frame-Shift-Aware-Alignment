"""
@author: Jules Kreuer, Catalina, Samira Breitling, Sebastian Meyer
"""

import argparse
import numpy as np
import blosum
import fasta


def DP():
    
    """
    Dynamic Programming 
    
    Input:
        dnaseq: fasta_object or string, DNA sequence
        aaseq:  fasta_object or string, Aminoacid sequence
        gap:    Int, Gap-penatly
        shift:  Int, Frameshift-penalty
        blosum: Dict, BLOSUM dictionary 

    Returns:
        socore: Int, Score of aligment
        alignment: String, 
            Alignment in following format:
            MIHPFISLV\RP
            MIHPFISLV-RP
            With \ denoting a backward-framshift
            With / denoting a forward-framshift 
    """

    blosum_handler = blosum.BLOSUM(62)

    x = "CATCATATCI"
    y = "HHIJ"

    len_x = len(x)
    len_y = len(y)

    alignment = np.zeros((len_y+3, len_x+1), dtype=int)
    alignment[0, :] = [x+2 for x in range(0, -(len_x+1), -1)]
    alignment[:, 0] = [y-2 for y in range(0, -(len_y+3), -1)]
    alignment[:, 1] = [y-1 for y in range(0, -(len_y + 3), -1)]
    alignment[:, 2] = [y for y in range(0, -(len_y + 3), -1)]

    print(alignment)

    for i in range (0, len_x):
        for j in range (0, len_y):

            translated_codon_i3 = fasta.translate_seq(x[i-2:i+1])
            translated_codon_i2 = fasta.translate_seq(x[i-1:i])
            translated_codon_i = fasta.translate_seq(x[i:i])
            translated_codon_i = fasta.translate_seq(x[i:i])

            #amino_acid = y[j]

            print(translated_codon, amino_acid)

            # Case MATCH
            if translated_codon == amino_acid:
                 = blosum_handler.get(translated_codon, amino_acid)

    return




def main():

    # TODO IO Functions
        # READ DNA
        # READ AA
        # READ BLOSUM

    # TODO MAIN DP
        # DNA Translation
        # AA Translation
        # Traceback
    DP()

    # TODO IO
        # PRINT Output
        # WRITE Output

    return

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description='A NICE TEXT')
    parser.add_argument('-d', "--dnaseq", action='store', dest='dnaseq', 
                        help='Specify path to DNA Sequence file (FASTA)', required=False)

    parser.add_argument('-aa', "--aaseq", action='store', dest='aaseq', 
                        help='Specify path to amino-acid file (FASTA)', required=False)

    parser.add_argument('-gp', "--gap", action='store', dest='gap', 
                        help='DNA gap penalty', required=False)

    parser.add_argument('-s', "--shift", action='store', dest='shift', 
                        help='DNA gap penalty', required=False)

    parser.add_argument('-b', "--blosum", action='store', dest='gap', 
                        help='Specify path to blosum matrix. Default: Blosum62', required=False)

    parser.add_argument('-o', "--out", action='store', dest='out', 
                        help='Specify path to output file.', required=False)

    parser.add_argument('-q', "--quiet", action='store_true', dest='quiet', 
                        help='Less verbose output', required=False)

    args = parser.parse_args()
    #FILE_PATH = args.path
    # TODO Read all args
    main()