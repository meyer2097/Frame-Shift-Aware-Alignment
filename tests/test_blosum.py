# Testing BLOSUM Handling
from frameshift_aware_alignment import BLOSUM
import pytest
from os import path


f = float("-inf")


@pytest.mark.parametrize("blosum_number,expected", [
    (45, [0, 1, -2, -5, 3, 1, -3, -5, -2, -2, 1, -5, f, f, f, f]),
    (50, [0, 1, -1, -5, 3, 2, -4, -5, -3, -1, 1, -5, f, f, f, f]),
    (62, [0, 0, -1, -4, 2, 1, -3, -4, -3, -2, 1, -4, f, f, f, f]),
    (90, [0, 1, -2, -6, 2, 1, -4, -6, -4, -3, 0, -6, f, f, f, f])
])
def test_blosum(blosum_number, expected):
    bm = BLOSUM(blosum_number)
    get_test = []

    for a in ["H", "K", "W", "U"]:
        for b in ["R", "Q", "F", "*"]:
            get_test.append(bm.get(a, b))

    assert get_test == expected


@pytest.mark.filterwarnings("ignore:Blosum")
def test_blosum_custom_file():
    fp = path.join(path.dirname(__file__), "test.blosum")
    bm = BLOSUM(fp)
    labels = ["A", "R", "N", "D"]
    s = sum([bm.get(a, b) for b in labels for a in labels])
    assert s == 17


@pytest.mark.xfail
def test_blosum_invalid_file():
    fp = path.join(path.dirname(__file__), "fail.blosum")
    bm = BLOSUM(fp)
    bm.get("A", "B")


@pytest.mark.xfail
def test_blosum_empty():
    bm = BLOSUM(None)
    bm.get("A", "B")
