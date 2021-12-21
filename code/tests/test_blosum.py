# Testing BLOSUM Handling
import pytest
from blosum import BLOSUM


def test_blosum_get():
    blosum_results = {
        45: [0, 1, -2, -5, 3, 1, -3, -5, -2, -2, 1, -5, float("-inf"), float("-inf"), float("-inf"), float("-inf")],
        50: [0, 1, -1, -5, 3, 2, -4, -5, -3, -1, 1, -5, float("-inf"), float("-inf"), float("-inf"), float("-inf")],
        62: [0, 0, -1, -4, 2, 1, -3, -4, -3, -2, 1, -4, float("-inf"), float("-inf"), float("-inf"), float("-inf")],
        90: [0, 1, -2, -6, 2, 1, -4, -6, -4, -3, 0, -6, float("-inf"), float("-inf"), float("-inf"), float("-inf")]
    }

    for idx, blosum_number in enumerate(blosum_results.keys()):
        blosum = BLOSUM(blosum_number)
        get_test = []

        for a in ["H", "K", "W", "U"]:
            for b in ["R", "Q", "F", "*"]:
                get_test.append(blosum.get(a, b))

        assert get_test == blosum_results[blosum_number]
