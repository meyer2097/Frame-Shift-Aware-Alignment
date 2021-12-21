# Testing BLOSUM Handling
from frameshift_aware_alignment import align


def test_optimal_align():
    dnaSeq = "CGGCCTTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    gap = 2
    shift = 15
    bm = 62
    score, dnaSeq_align, aaSeq_align = align(dnaSeq, aaSeq, gap, shift, bm)
    assert aaSeq == dnaSeq_align
    assert aaSeq == aaSeq_align
    assert score == 32


def test_forward_align():
    dnaSeq = "CGG" + "T" + "CCTTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    gap = 2
    shift = 15
    bm = 62
    score, dnaSeq_align, aaSeq_align = align(dnaSeq, aaSeq, gap, shift, bm)
    assert "R/PSIFF" == dnaSeq_align
    assert "R-PSIFF" == aaSeq_align
    assert score == 17


def test_backwards_align():
    dnaSeq = "CGG" + "TT" + "CCTTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    gap = 2
    shift = 15
    bm = 62
    score, dnaSeq_align, aaSeq_align = align(dnaSeq, aaSeq, gap, shift, bm)
    assert "R\\PSIFF" == dnaSeq_align
    assert "R-PSIFF" == aaSeq_align
    assert score == 17


def test_incodon_align():
    dnaSeq = "CGGCC" + "G" + "TTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    gap = 2
    shift = 15
    bm = 62
    score, dnaSeq_align, aaSeq_align = align(dnaSeq, aaSeq, gap, shift, bm)
    assert "RP/SIFF" == dnaSeq_align
    assert "RP-SIFF" == aaSeq_align
    assert score == 17
