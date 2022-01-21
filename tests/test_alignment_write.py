# Testing writing of Alignment
import frameshift_aware_alignment as faa
from os import remove

gap = 2
shift = 15
bm = 62


def test_write():
    dnaSeq = "CGG" + "T" + "CCTTCTATCTTCTTC"
    aaSeq = "RPSIFF"
    score, dnaSeq_align, aaSeq_align = faa.align(dnaSeq, aaSeq, gap, shift, bm, out="test_write")
    assert "R/PSIFF" == dnaSeq_align
    assert "R-PSIFF" == aaSeq_align
    assert score == 17

    with open("test_write", "r") as f:
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]
        assert "DNA: " + dnaSeq_align in lines
        assert "AA : " + aaSeq_align in lines

    remove("test_write")
