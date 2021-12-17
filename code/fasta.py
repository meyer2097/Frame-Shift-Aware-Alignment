"""
@author: Jules Kreuer
A simple fasta object, reader and writer 
"""

from os import path

# Translation Dictionary according to 
# https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
translation_dict = {"TTT": "F", "TTC": "F",
                    "TTA": "L", "TTG": "L",
                    "TCT": "S", "TCC": "S","TCA": "S", "TCG": "S",
                    "TAT": "Y", "TAC": "Y",
                    "TAA": "*", "TAG": "*", "TGA": "*",
                    "TGT": "C", "TGC": "C",
                    "TGG": "W", 
                    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
                    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                    "CAT": "H", "CAC": "H",
                    "CAA": "Q", "CAG": "Q",
                    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                    "ATT": "I", "ATC": "I", "ATA": "I",
                    "ATG": "M", 
                    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                    "AAT": "N", "AAC": "N",
                    "AAA": "K", "AAG": "K",
                    "AGT": "S", "AGC": "S",
                    "AGA": "R", "AGG": "R",
                    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
                    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                    "GAT": "D", "GAC": "D",
                    "GAA": "E", "GAG": "E",
                    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

complement_dict = { "A": "T",
                    "G": "C",
                    "C": "G",
                    "T": "A",
                    "U": "A",
                    "R": "Y",
                    "Y": "R",
                    "S": "S",
                    "W": "W",
                    "K": "M",
                    "M": "K",
                    "B": "V", 
                    "V": "B", 
                    "D": "H", 
                    "H": "D"}


class fasta_object():
    def __init__(self, head, body):
        """
        Object to keep a valid fasta object
        """
        if head.startswith(">"):
            self.head = head
        else:
            self.head = f">{head}"
        
        self.body = body

    def __str__(self): 
        """
        Magic method to allow fasta_object printing.
        """
        out = f'{self.head}\n'

        # Print only 70 chars per line
        for i in range(0, len(self.body), 70):
            out += f'{self.body[i:i+70]}\n'
        return out[:-1] # Remove tailing newline

    def __repr__(self):
        """
        Magic method to allow printing of fasta_object representation.
        """
        return f'fasta_object("{self.head}", "{self.body}")'

    def __eq__(self, o):
        """
        Magic method to allow equality check on fasta_objects.
        Does not check for header equality.
        """
        return self.body == o.body
    
    def toAmino(self, d=translation_dict):
        """
        Translates the dna sequence of a fasta_object to amino-acids.
        Reading frame starts at position 0, tailing bases will be ignored.
        Attention: Will throw exception if triplet is not found.
        """
        self.body = translate_seq(self.body, d)
    
    def toReverseComp(self):
        """
        Translates the dna sequence the reverse complement.
        """
        self.body = reverse_comp(self.body)

def read_fasta(fileName):
    """
    Reads a fasta-style file and returns a list of fasta_objects
    """

    if not path.isfile(fileName):
        raise FileNotFoundError("Fasta File not found!")

    fasta_objects = []
    with open(fileName, 'r') as f:
        head = ""
        body = ""
        newObject = True
        for line in f:
            # go through each line
            # First Header
            if newObject and line.startswith(">"):
                head = line.strip()
                body = ""
                newObject = False
            # N-th Header
            elif line.startswith(">"):
                fasta_objects.append(fasta_object(head, body))
                head = line.strip()
                body = ""
            # Sequence
            else:
                body += line.strip().upper()
        # append last element
        fasta_objects.append(fasta_object(head, body))
    return fasta_objects

def write_fasta(fasta_pairs, fileName, mode="w"):
    """
    Writes a list of fasta_objects or a single one to a file.
    Takes fasta_objects as input.
    """

    if not isinstance(fasta_pairs, list):
        fasta_pairs = [fasta_pairs]

    with open(fileName, mode) as f:
        for fo in fasta_pairs:
            f.write(f"{fo.head}\n")
            body_len = len(fo.body)
            # Write only 70 chars per line
            for i in range(0, body_len, 70):
                f.write(f"{fo.body[i:i+70]}\n")
    return True

def print_fasta(fasta):
    """
    Prints a single or a list of fasta_objects.
    """

    if not isinstance(fasta, list):
        fasta = [fasta]
    
    for fo in fasta:
        print(fo.head)
        body_len = len(fo.body)
        # Print only 70 chars per line
        for i in range(0, body_len, 70):
            print(fo.body[i:i+70])
    return None

def __maybeFind(key, d, alt):
    try:
        return d[key]
    except KeyError:
        return alt

def translate_seq(seq, d=translation_dict):
    """
    Translates a DNA sequence to a AA sequence.
    Reading frame starts at position 0, tailing bases will be ignored.
    Attention: Will throw exception if triplet is not found.

    To translate a fasta_object use objcet.toAmino()
    
    Input:
        seq: String, sequence to translate
        d: dict, dictionary of translation
    Returns:
        translated: String, translated sequence
    """

    translated = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if not len(codon) == 3:
            break
        translated += __maybeFind(codon, d, "~")
    return translated



def reverse_comp(seq):
    """
    Reverses complement of sequence.
    """
    seq = reversed(seq)
    seq = "".join(map(lambda b: __maybeFind(b, complement_dict, b), seq))
    return seq