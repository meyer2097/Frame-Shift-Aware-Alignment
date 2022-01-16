# Testing BLOSUM Handling
import frameshift_aware_alignment as faa

gap = 2
shift = 15
bm = 62


def test_optimal_align():
    dnaSeq = "CGGCCTTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    score, dnaSeq_align, aaSeq_align = faa.align(dnaSeq, aaSeq, gap, shift, bm)
    assert aaSeq == dnaSeq_align
    assert aaSeq == aaSeq_align
    assert score == 32


def test_forward_align():
    dnaSeq = "CGG" + "T" + "CCTTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    score, dnaSeq_align, aaSeq_align = faa.align(dnaSeq, aaSeq, gap, shift, bm)
    assert "R/PSIFF" == dnaSeq_align
    assert "R-PSIFF" == aaSeq_align
    assert score == 17


def test_backwards_align():
    dnaSeq = "CGG" + "TT" + "CCTTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    score, dnaSeq_align, aaSeq_align = faa.align(dnaSeq, aaSeq, gap, shift, bm)
    assert "R\\PSIFF" == dnaSeq_align
    assert "R-PSIFF" == aaSeq_align
    assert score == 17


def test_incodon_align():
    dnaSeq = "CGGCC" + "G" + "TTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    score, dnaSeq_align, aaSeq_align = faa.align(dnaSeq, aaSeq, gap, shift, bm)
    assert "RP/SIFF" == dnaSeq_align
    assert "RP-SIFF" == aaSeq_align
    assert score == 17


def test_deletion_dna_align():
    dnaSeq = "CGGTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    score, dnaSeq_align, aaSeq_align = faa.align(dnaSeq, aaSeq, gap, shift, bm)
    assert "R-SIFF" == dnaSeq_align
    assert "RPSIFF" == aaSeq_align
    assert score == 23


def test_deletion_aa_align():
    dnaSeq = "CGGCCTTCTATCTTCTTC"
    aaSeq = "RSIFF"
    score, dnaSeq_align, aaSeq_align = faa.align(dnaSeq, aaSeq, gap, shift, bm)
    assert "RPSIFF" == dnaSeq_align
    assert "R-SIFF" == aaSeq_align
    assert score == 23
