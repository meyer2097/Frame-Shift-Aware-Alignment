import argparse
from src.frameshift_aware_alignment.frameshift_aware_alignment import align


if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="""
A frameshift aware Needleman-Wunsch-Global-Aligner for DNA- and amino acid-sequences.
""")
    parser.add_argument('-d', "--dnaseq", action='store', dest='dnaseq',
                        help="DNA sequence or path to a single entry fasta file.",
                        required=True)

    parser.add_argument('-aa', "--aaseq", action='store', dest='aaseq',
                        help="Aminoacid sequence or path to a single entry fasta file.",
                        required=True)

    parser.add_argument('-gp', "--gap", action='store', dest='gap',
                        help="Gap penalty (absolute value)",
                        required=True)

    parser.add_argument('-sp', "--shift", action='store', dest='shift',
                        help="Frameshift penalty (absolute value)",
                        required=True)

    parser.add_argument('-b', "--blosum", action='store', dest='blosum', default=62,
                        help=""" Specify blosum matrix. One of: 45, 50, 62, 80, 90. Default: 62.
                                 Will be irgnored if -bp is set.""",
                        required=False)

    parser.add_argument('-bp', "--blosum_path", action='store', dest='blosum_path',
                        help="""Specify path to a custom blosum matrix.
Use only if -b is not a viable solution.""",
                        required=False)

    parser.add_argument('-o', "--out", action='store', dest='out', default=False,
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
