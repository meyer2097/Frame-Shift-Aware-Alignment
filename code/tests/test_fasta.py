import pytest
import fasta
from os import path, remove

def test_read_fasta():
    fp = path.join(path.dirname(__file__), "read_fasta_test.fasta")
    assert fasta.read_fasta(fp) == [
        fasta.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
        fasta.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]


def test_write_read_fasta():

    file_path = path.join(path.dirname(__file__), "write_fasta_test.fasta")

    fo = [fasta.fasta_object(">Atlantic dolphin", "CGGCCTT*CTAAAAATTZZZ*ZZZZASASD*TCTTCTTC"),
          fasta.fasta_object(">Pacific dolphin", "CTTTCTATCTCSATTTCCTCT")]
    
    fasta.write_fasta(fo, file_path)
    
    fo_read = fasta.read_fasta(file_path)
    remove(file_path)
    assert fo == fo_read
    


def test_translate_seq():
    assert fasta.translate_seq("CGGCCTTCTATCTTCTTC") == "RPSIFF"
    assert fasta.translate_seq("HELLO") == "incorrect Sequence"


def test_print_fasta(capsys):
    file_path = path.join(path.dirname(__file__), "read_fasta_test.fasta")

    assert fasta.print_fasta(fasta.read_fasta(file_path)) is None

    # check if it prints the sequences correctly
    captured = capsys.readouterr()
    assert captured.out == ">Atlantic dolphin" + "\n" + "CGGCCTTCTATCTTCTTC" + "\n" + ">Pacific dolphin" + "\n" + "CTTTCTATCTCTTTCCTCT" + "\n"
    assert fasta.print_fasta(fasta.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC")) is None
    # prints: ">Atlantic dolphin" + \"\n" + "CGGCCTTCTATCTTCTTC"
    assert fasta.print_fasta([fasta.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
                              fasta.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]) is None
    # prints ">Atlantic dolphin" + "\n" + "CGGCCTTCTATCTTCTTC" + "\n" + ">Pacific dolphin" + "\n" +
    # "CTTTCTATCTCTTTCCTCT"
