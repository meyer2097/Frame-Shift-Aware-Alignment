"""
Frameshift aware Needleman-Wunsch for DNA and AA-Sequences.
"""
from .frameshift_aware_alignment import align
from .fasta import fasta_object, read_fasta, write_fasta, print_fasta, translate_seq
from .blosum import BLOSUM
__all__ = [
    "align",
    "fasta_object",
    "read_fasta",
    "write_fasta",
    "print_fasta",
    "translate_seq",
    "BLOSUM"
]
