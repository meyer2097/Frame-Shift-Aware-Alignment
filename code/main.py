"""
@author: Jules Kreuer, Catalina, Samira Breitling, Sebastian Meyer
"""

import argparse
import numpy as np
import blosum
import fasta


def DP(shift=-1, gap=-1):
    
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

    x = "CATCATATC"
    y = "HHI"

    len_x = len(x)
    len_y = len(y)

    alignment = np.zeros((len_x+3, len_y+1), dtype=int)
    alignment[:, 0] = [x+2 for x in range(0, -(len_x+3), -1)]
    alignment[0, :] = [y-2 for y in range(0, -(len_y+1), -1)]
    alignment[1, :] = [y-1 for y in range(0, -(len_y + 1), -1)]
    alignment[2, :] = [y for y in range(0, -(len_y + 1), -1)]
    alignment[0,0]=0
    alignment[0,1]=0
    alignment[1,0]=0
    alignment[1,1]=0
    #print(alignment)
    
    for i in range (3, len_x+3):
        for j in range (1, len_y+1):
            
            
            translated_codon_i = fasta.translate_seq(x[i-2:i+1])
            translated_codon_i1 = fasta.translate_seq(x[i-3:i])
            translated_codon_i2 = fasta.translate_seq(x[i-4:i-1])
            translated_codon_i3 = fasta.translate_seq(x[i-5:i-2])
            
            acid_j1= y[j-1]
            acid_j= y[-1]
            
            case1=alignment[i-3,j-1]+blosum_handler.get(translated_codon_i, acid_j)
            case2=alignment[i,j-1]+blosum_handler.get(translated_codon_i, acid_j1)+shift
            case3=alignment[i-1,j]+blosum_handler.get(translated_codon_i1, acid_j)+gap
            case4=alignment[i-2,j]+blosum_handler.get(translated_codon_i2, acid_j)+2*gap
            case5=alignment[i-3,j]+blosum_handler.get(translated_codon_i3, acid_j)+shift
            
            alignment[i,j]=max(case1,case2,case3,case4,case5)
            #amino_acid = y[j]
    print(alignment)
           

            

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