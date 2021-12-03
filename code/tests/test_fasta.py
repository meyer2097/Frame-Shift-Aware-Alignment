import pytest
import fasta


def test_read_fasta():
    assert fasta.read_fasta("read_fasta_test.fasta") == [
        fasta.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
        fasta.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]


def test_write_fasta():
    fasta.write_fasta([fasta.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
                       fasta.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")], "write_fasta_test.fasta")
    file_test = open("write_fasta_test.fasta", "r")
    file = open("read_fasta_test.fasta", "r")
    assert file_test.read() == file.read() + "\n"


def test_translate_seq():
    assert fasta.translate_seq("CGGCCTTCTATCTTCTTC") == "RPSIFF"
    assert fasta.translate_seq("HELLO") == "incorrect Sequence"


def test_print_fasta(capsys):
    assert fasta.print_fasta(fasta.read_fasta("read_fasta_test.fasta")) is None
    # check if it prints the sequences correctly
    captured = capsys.readouterr()
    assert captured.out == ">Atlantic dolphin" + "\n" + "CGGCCTTCTATCTTCTTC" + "\n" + ">Pacific dolphin" + "\n" + "CTTTCTATCTCTTTCCTCT" + "\n"
    assert fasta.print_fasta(fasta.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC")) is None
    # prints: ">Atlantic dolphin" + \"\n" + "CGGCCTTCTATCTTCTTC"
    assert fasta.print_fasta([fasta.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
                              fasta.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]) is None
    # prints ">Atlantic dolphin" + "\n" + "CGGCCTTCTATCTTCTTC" + "\n" + ">Pacific dolphin" + "\n" +
    # "CTTTCTATCTCTTTCCTCT"
