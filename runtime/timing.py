# Measure runtime
from frameshift_aware_alignment import __alignment_core, BLOSUM
import random
from timeit import repeat
from statistics import mean, median
import pandas as pd

gap = 2
shift = 15
bm = 62

alphaD = "AGTC"
alphaA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
          "M", "F", "P", "S", "T", "W", "Y", "V", "B", "J", "Z", "X"]


def rSeq(alpha, n):
    return ''.join(random.choice(alpha) for x in range(n))


def a(i, j):
    __alignment_core(rSeq(alphaD, i), rSeq(alphaA, j), 5, 10, bm)


bigN = 20
dna_low = 1
aa_low = 1
dna_high = 180
aa_high = 60
bm = BLOSUM(62)


if __name__ == "__main__":
    timeDict = {}
    for i in range(dna_low, dna_high+1):
        print(f"Currently at: {i}")
        timeDict[i] = {}
        for j in range(aa_low, aa_high+1):
            timeDict[i][j] = repeat(lambda: a(i, j), repeat=bigN, number=1)

    # Write median to csv file
    medianDict = {}
    for i in timeDict:
        medianDict[i] = {}
        for j in timeDict[i]:
            medianDict[i][j] = median(timeDict[i][j])

    df = pd.DataFrame.from_dict(medianDict, orient='index')
    with open("median.csv", "w") as f:
        lines = df.to_csv(index=True)
        f.writelines(lines)

    meanDict = {}
    for i in timeDict:
        meanDict[i] = {}
        for j in timeDict[i]:
            meanDict[i][j] = mean(timeDict[i][j])

    df = pd.DataFrame.from_dict(meanDict, orient='index')
    with open("mean.csv", "w") as f:
        lines = df.to_csv(index=True)
        f.writelines(lines)
