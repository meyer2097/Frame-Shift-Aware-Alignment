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
In an alignment the gap penalty is an important parameter to score alignments. In our frame-shift aware alignment program we first considered the linear gap penalty, which penalizes each gap equally. However, this is not realistic when it comes to biological events. Each mutation (insertion or deletion) is an indivitiual evolutionary event. Consequently we consider an affine gap model which penalizes a gap opening higher than a gap extend. This gap model is called affine gap penalty. 

Linear gap penalty: AAGTGTGCGTTCCGATT
                    AA--GT--G---CGATT
                    
Affine gap penalty: AAGTGTGCGTTCCGATT
                    AAGTG-------CGATT
                    
In addition to the plain frame-shift aware alignment, we also implemented a frame-shift alignment which also considers the affine gap penalty model.


## Usage
Usage as module:

```
# Clone repo
cd Frame-Shift-Aware-Alignment
pip install .
```

```python
from frameshift-aware-alignment import align

score, dna, aa = align(...)
```

Usage as standalone script
```
python3 align.py [-h] -d DNASEQ -aa AASEQ -gp GAP -sp SHIFT [-b BLOSUM] [-bp BLOSUM_PATH] [-o OUT]
```


