import numpy as np

def __alignment_core(dnaSeq: str, aaSeq: str,
                     gap: int, gap_open:int, gap_extend:int, frameshift:int):
    """
    Core implementation of the frameshift aware Needleman-Wunsch for DNA and AA-Sequences.

    Parameters:
        dnaSeq:  string, dna sequence.
        aaSeq:   string, amino-acid sequence.
        gap:     int, score for gap.
        shift:   int, score for frameshift.
        bm:      BLOSUM Object, blosum matrix for AA-alignment scores.
    Returns:
        score: Int, Score of aligment
        dnaSeq_align: String, alignment of the DNA sequence.
            With \\ denoting a backward-framshift.
            With / denoting a forward-framshift.
            With - denoting a gap.
        aaSeq_align: String, alignment of the AA sequence.
            With - denoting a gap.
    """

    n = len(aaSeq) + 1  # Row
    m = len(dnaSeq) + 3  # Col

    # Add space to compensate index shifting for amino-acid
    aaSeq = f" {aaSeq}"

    # Init matrix with zeros
    score_matrix = np.zeros((n, m), dtype=float)
    insertion_matrix = np.zeros((n, m), dtype=float)
    deletion_matrix = np.zeros((n, m), dtype=float)


    # Init traceback matrix
    traceback_matrix = np.zeros((n, m), dtype=str)

    # First three columns
    """for i in range(n):
        score_matrix[i][0] = -i*gap
        traceback_matrix[i][0] = "U"

        score_matrix[i][1] = -i*gap - shift
        traceback_matrix[i][1] = "1"

        score_matrix[i][2] = -i*gap - shift
        traceback_matrix[i][2] = "2"""

    score_matrix[0][0] = float("-inf")
    deletion_matrix[0][0] = float("-inf")
    insertion_matrix[0][0] = float("-inf")
    score_matrix[1][0] = float("-inf")
    deletion_matrix[1][0] = float("-inf")

    g = n
    print(g)

    for g_run in range(1, g+1):
        deletion_matrix[3*g_run + 1][0] = -gap_open-gap_extend*g_run
        """insertion_matrix[1][g_run] = -gap_open-gap_extend*g_run
        deletion_matrix[3 * g_run + 1][0] = -gap_open - gap_extend * g_run
        deletion_matrix[3 * g_run][0] = -gap_open - gap_extend * g_run - gap
        deletion_matrix[3 * g_run + 2][0] = -gap_open - gap_extend * g_run - gap

        score_matrix[0][g_run] = float("-inf")
        score_matrix[1][g_run] = float("-inf")
        score_matrix[2][g_run] = float("-inf")
        score_matrix[3][g_run] = float("-inf")

        deletion_matrix[0][g_run] = float("-inf")
        deletion_matrix[1][g_run] = float("-inf")
        deletion_matrix[2][g_run] = float("-inf")
        deletion_matrix[3][g_run] = float("-inf")

        score_matrix[g_run + 1][0] = float("-inf")

        insertion_matrix[0][g_run] = float("-inf")
        insertion_matrix[g_run + 1][0] = float("-inf")"""

    print(score_matrix)
    print(deletion_matrix)
    print(insertion_matrix)

    for i in range(3, m):
        for j in range(1, n):
            translated_seq = ft.translate_seq(dnaSeq[j - 3:j])

            s_max_1 = max(score_matrix[i-3][j-1], deletion_matrix[i-3][j-1], insertion_matrix[i-3][j-1]) + translated_seq
            s_max_2 = max(score_matrix[i-3][j-1], deletion_matrix[i-3][j-1], insertion_matrix[i-3][j-1], score_matrix[i-4][j-1], deletion_matrix[i-4][j-1], insertion_matrix[i-4][j-1]) - frameshift + translated_seq

            s_max = max(s_max_1, s_max_2)

            score_matrix[i][j] = s_max

            i_max = max(score_matrix[i][j-1]-gap_open, insertion_matrix[i][j-1]) - gap_extend
            insertion_matrix[i][j] = i_max


            d_max_1 = max(score_matrix[i-3][j] - gap_open, deletion_matrix[i-3][j]) - gap_extend
            d_max_2 = max(score_matrix[i-2][j] - gap_open, deletion_matrix[i-2][j], score_matrix[i-4][j] - gap_open, deletion_matrix[i-4][j]) - frameshift - gap_extend
            d_max = max(d_max_1, d_max_2)

            deletion_matrix[i][j] = d_max


__alignment_core("GGGGGG", "VV", 3, 5, 10)