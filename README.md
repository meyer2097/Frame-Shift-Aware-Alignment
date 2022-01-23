# Frame-Shift-Aware-Alignment
Sequence Bioinformatics Group Project WS 21/22
Samira Breitling, Jules Kreuer, Sebastian Meyer and Catalina Schlotterer

## Introduction
Sequence alignment is the basis of many bioinformatic problems. The Needleman-Wunsch alorithm as well as the Smith-Waterman algorithm are dynamic programming approaches to compare two sequences and compute a global/local aligment. These algorithm perform reliable alignments for short aas well as high quality reads. However, when comparing a DNA sequence to an amino acid sequence erronous insertions and deletions result in so called frameshifts. The following example illustrates this issue:

When inserting three bases in a error-free read:

![inserting bases to DNA strand lead to frameshift](https://user-images.githubusercontent.com/94982104/150514625-b0d299fa-2a21-4fc0-b1f8-5c53d9637f06.png)

Then translate this read to its aminoacid sequence and compare it with the actual protein:

![inserting bases to DNA strand lead to frameshift](https://user-images.githubusercontent.com/94982104/150515030-3abdf8ed-db21-4f92-ab7a-332c8111fabc.png)

The green part of the alignment show perfect matching parts of the alignment. However between the first and third insertion, where the alignment is again in the original frame, the alignment is of poor quality.
A frameshift aware alignment deals with this frameshifts and reports forward frameshifts (inserting two bases) and backward frameshifts (inserting one basis) in the alignment, which leads to a alignment of significantly higher alignment quality.

![inserting bases to DNA strand lead to frameshift](https://user-images.githubusercontent.com/94982104/150515892-70ba849d-47aa-4007-9f55-c0ec0efbbb3d.png)

In this group project we constructed a dynamic programm which reports an alignment and the corresponding score considering frameshifts.

### Affine gap penalty
In an alignment the gap penalty is an important parameter to score alignments. In our frameshift aware alignment program we first considered the linear gap penalty, which penalizes each gap equally. However, this is not realistic when it comes to biological events. Each mutation (insertion or deletion) is an indivitiual evolutionary event. Consequently we consider an affine gap model which penalizes a gap opening higher than a gap extend. This gap model is called affine gap penalty. 

![affine-gap-penalty-vs-linear](https://user-images.githubusercontent.com/94982104/150522224-2640f306-3508-4014-b41b-c03c2539e0da.png)
                    
In addition to the plain frameshift aware alignment, we also implemented a frameshift alignment which also considers the affine gap penalty model.


## Teamwork
In the beginning of the group project we first implemented the different helper functions: translating the DNA sequence to an amino acid sequence, reading in the fasta files and reading in the blosum matrix. We tested these functions and added different blossum matrices. Additionally we added the main function and discussed the input paramteters and the output format. Next we started with the actual frameshift aware alginment. Therefore we decided to implement two algorithms in pairs of two and compared both implementations. This allowed us to have two different approaches where we could discuss difficulties and pair programming enabled efficient implementing of the algorithm. After the implementation of the frameshift aware alignment, we considered integrating the affine gap penalty. After some research we came up with a paper that seemed to solve this problem. However, they used a different approach doing the frameshift alignment than we did. Consequently, one group decided to implement the affine gap penalty by adding the affine gap penalty restrictions to our affine gap penalty implementation. The other group tried to implement as described in the paper. However, this approach seemed to be too challenging which is why we decided for the other implementation. Therefore working in pairs to resolve the tasks was a good idea.


## Algorithm
As already mentioned in the introduction the frame-shift aware algorithm is a dynamic porgramm. Hence, a matrix which scores is filled to compute the alignment and the score. In the picture below it is illustrated how to compute the score for a a cell. The black errors indicate the same cases as in a "basic" alignment, while the blue errors show which cell we need to consider for frameshifts. In this case we either insert one base or two bases into the DNA query sequence q.

<img src="https://user-images.githubusercontent.com/94982104/150526888-85e380a5-543e-48be-998a-0cef56917dba.png" style="width:700px; max-width:90%" alt="table img">

<img src="https://user-images.githubusercontent.com/94982104/150528688-b96e248d-fc5d-405b-a293-e2b56e740526.png" style="width:600px; max-width:90%" alt="rescursion img">
where d is the gap penalty, and k the frameshift penalty. 

Note that the algorithm only works if the frameshift penalty is larger than the gap penalty. Hence, the occurance of an insertion or deletion is more probable than the occurance of a frameshift.

### Affine gap penalty
Several papers have shown that an implementation of a frameshift awaren aligner with affine gap penalty is possible. For this project, attempts were made to implement this. However, it quickly became apparent that this was too big an undertaking. 

The following section shows the **incomplete** attempt. The code can be found in the branch 'affine_test_jules'.

For the affine gap penalty we need two additional matrices:
The top matrix should provide a minimal score if the extension of an amino-acid insertion is optimal, whereas the bottom matrix is responsible for the dna insertions.

We initialized the three matrices as follows:

_Iniitialization columns_
```python
basePenalty = -3*gap_open - (i*gap_extend)

score_matrix[i][0] = basePenalty
top_matrix[i][0] = top_matrix[i][1] = top_matrix[i][2] = - infinity
score_matrix[i][1] = score_matrix[i][2] = basePenalty - k
```
_Initialization rows_
```python
basePenalty = -3*gap_open - (j*gap_extend)

score_matrix[0][j] = basePenalty
score_matrix[0][j+1] = score_matrix[0][j+2] = basePenalty - k
bottom_matrix[0][j] = bottom_matrix[0][j+1] = bottom_matrix[0][j+2] = - infinity
```
Recursion:
<img src="https://user-images.githubusercontent.com/94982104/150538090-9711a9e7-48b2-4219-ad37-9d78b5ed75fa.png" alt="rescursion 2 img">


_All figures except for the affine gap penalty recursion were taken out of the Sequence Bioinformatics script by Daniel Huson_

## Runtime
![Runtime](https://user-images.githubusercontent.com/25013642/150572224-4b660955-0101-45ac-8c61-cede61ae3280.png)
The frameshift-aware global alginer has a theretical runtime `O(AA*DNA)` with AA=length of the amino-acid sequence and DNA=length of dna sequence.

To determine the actual runtime, random DNA and AA sequences with different length were genereated and aligned. In order to obtain more accurate results, each length combination was repeated 20 times, the median was taken.
The results can be viewed in the previous heatmap.

To predict the running time, the following formula, determined by a simple linear model, can be used:

`time_in_ms = 0.70030  + 0.00778 * AA * DNA`

## Installation
This project was build a python module and can be installed with pip:

```
# Clone repo
cd Frame-Shift-Aware-Alignment
pip install .
```

### Usage
To align two sequences, the function `align` can be called. Some sequences for testing can be found within the `tests` folder.
Several arguments should be provided:
```python
    Parameters:
        dnaSeq:  string, dna sequence.
        aaSeq:   string, amino-acid sequence.
        gap:     int, score for gap.
        shift:   int, score for frameshift.
        bm:      int, string, One of [45, 50, 62, 80, 90] or path to blosum matrix.
        verbose: bool, print score and alignment.
    Returns:
        score: Int, Score of aligment
        dnaSeq_align: String, alignment of the DNA sequence.
            With \\ denoting a backward-framshift.
            With / denoting a forward-framshift.
            With - denoting a gap.
        aaSeq_align: String, alignment of the AA sequence.
            With - denoting a gap.
```
**Example**
```python
from frameshift_aware_alignment import align

score, dna, aa = align("AGTCAGT", "MM", 6, 15, 62, false)
```

## Standalone version
A usage without installation is also possible. After downloading the arcive, execute top-level file `align.py`.

```
python3 align.py [-h] -d DNASEQ -aa AASEQ -gp GAP-Penalty -sp SHIFT-Penalty [-b BLOSUM] [-bp BLOSUM_PATH] [-o OUT]
```


