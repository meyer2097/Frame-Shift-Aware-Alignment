"""
@author: Jules Kreuer
"""
import argparse
import fasta_reader

# Standard translation table according to https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
translation_dict = {"GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
                    "AAC": "B", "AAT": "B", "GAC": "B", "GAT": "B",
                    "TGC": "C", "TGT": "C",
                    "GAC": "D", "GAT": "D",
                    "GAA": "E", "GAG": "E",
                    "TTC": "F", "TTT": "F",
                    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
                    "CAC": "H", "CAT": "H",
                    "ATA": "I", "ATC": "I", "ATT": "I",
                    "AAA": "K", "AAG": "K",
                    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L", "TTA": "L", "TTG": "L",
                    "ATG": "M",
                    "AAC": "N", "AAT": "N",
                    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
                    "CAA": "Q", "CAG": "Q",
                    "AGA": "R", "AGG": "R", "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
                    "AGC": "S", "AGT": "S", "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
                    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
                    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
                    "TGG": "W",
                    "NNN": "X",
                    "TAC": "Y", "TAT": "Y",
                    "CAA": "Z", "CAG": "Z", "GAA": "Z", "GAG": "Z",
                    "TAA": "*", "TAG": "*", "TGA": "*"}

        
def translate(fo):
    """
    Translate a single fasta_object.
    """
    translated = ""
    for i in range(0, len(fo.body), 3):
        codon = fo.body[i:i+3]
        if not len(codon) == 3:
            break 
        translated += translation_dict[codon]

    return fasta_object(fo.head + " translated", translated)
