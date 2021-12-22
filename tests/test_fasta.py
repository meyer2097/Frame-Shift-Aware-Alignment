# Test fasta reading, writing and translating
from frameshift_aware_alignment import fasta as ft
from os import path, remove


def test_read_fasta():
    fp = path.join(path.dirname(__file__), "read_fasta_test.fasta")
    assert ft.read_fasta(fp) == [
        ft.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
        ft.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]


def test_write_read_fasta():

    file_path = path.join(path.dirname(__file__), "write_fasta_test.fasta")

    fo = [ft.fasta_object(">Atlantic dolphin", "CGGCCTT*CTAAAAATTZZZ*ZZZZASASD*TCTTCTTC"),
          ft.fasta_object(">Pacific dolphin", "CTTTCTATCTCSATTTCCTCT")]

    ft.write_fasta(fo, file_path)

    fo_read = ft.read_fasta(file_path)
    remove(file_path)
    assert fo == fo_read


def test_translate_seq():
    assert ft.translate_seq("CGGCCTTCTATCTTCTTC") == "RPSIFF"
    assert ft.translate_seq("HELLO") == "~"


def test_print_fasta(capsys):
    file_path = path.join(path.dirname(__file__), "read_fasta_test.fasta")

    assert ft.print_fasta(ft.read_fasta(file_path)) is None

    # check if it prints the sequences correctly
    captured = capsys.readouterr()
    assert captured.out == """>Atlantic dolphin
CGGCCTTCTATCTTCTTC
>Pacific dolphin
CTTTCTATCTCTTTCCTCT
"""
    assert ft.print_fasta(ft.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC")) is None
    assert ft.print_fasta([ft.fasta_object(">Atlantic dolphin", "CGGCCTTCTATCTTCTTC"),
                           ft.fasta_object(">Pacific dolphin", "CTTTCTATCTCTTTCCTCT")]) is None
