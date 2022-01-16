# Test fasta reading, writing and translating
import frameshift_aware_alignment as faa
from os import path, remove


def test_read_fasta():
    fp = path.join(path.dirname(__file__), "read_fasta_test.fasta")
    assert faa.read_fasta(fp) == [
        faa.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
        faa.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]


def test_write_read_fasta():

    file_path = path.join(path.dirname(__file__), "write_fasta_test.fasta")

    fo = [faa.fasta_object(">Atlantic dolphin", "CGGCCTT*CTAAAAATTZZZ*ZZZZASASD*TCTTCTTC"),
          faa.fasta_object(">Pacific dolphin", "CTTTCTATCTCSATTTCCTCT")]

    faa.write_fasta(fo, file_path)

    fo_read = faa.read_fasta(file_path)
    remove(file_path)
    assert fo == fo_read


def test_translate_seq():
    assert faa.translate_seq("CGGCCTTCTATCTTCTTC") == "RPSIFF"
    assert faa.translate_seq("HELLO") == "~"


def test_print_fasta(capsys):
    file_path = path.join(path.dirname(__file__), "read_fasta_test.fasta")

    assert faa.print_fasta(faa.read_fasta(file_path)) is None

    # check if it prints the sequences correctly
    captured = capsys.readouterr()
    assert captured.out == """>Atlantic dolphin
CGGCCTTCTATCTTCTTC
>Pacific dolphin
CTTTCTATCTCTTTCCTCT
"""
    assert faa.print_fasta(faa.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC")) is None
    assert faa.print_fasta([faa.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
                           faa.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]) is None
