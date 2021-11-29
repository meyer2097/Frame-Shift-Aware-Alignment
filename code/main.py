"""
@author: Jules Kreuer, Catalina, Samira Breitling, Sebastian
"""

import argparse
import blosum



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

    # TODO IO
        # PRINT Output
        # WRITE Output

    return

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description='A NICE TEXT')
    parser.add_argument('-d', "--dnaseq", action='store', dest='dnaseq', 
                        help='Specify path to DNA Sequence file (FASTA)', required=True)

    parser.add_argument('-aa', "--aaseq", action='store', dest='aaseq', 
                        help='Specify path to amino-acid file (FASTA)', required=True)

    parser.add_argument('-gp', "--gap", action='store', dest='gap', 
                        help='DNA gap penalty', required=True)

    parser.add_argument('-s', "--shift", action='store', dest='shift', 
                        help='DNA gap penalty', required=False)

    parser.add_argument('-b', "--blosum", action='store', dest='gap', 
                        help='Specify path to blosum matrix. Default: Blosum62', required=False)

    parser.add_argument('-o', "--out", action='store', dest='out', 
                        help='Specify path to output file.', required=False)

    parser.add_argument('-q', "--quiet", action='store_true', dest='quiet', 
                        help='Less verbose output', required=False)

    args = parser.parse_args()
    FILE_PATH = args.path
    # TODO Read all args
    main()