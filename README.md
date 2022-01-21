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


_The figures were taken out of the Sequence Bioinformatics script by Daniel Huson_

## Teamwork
In the beginning of the group project we first implemented the different helper functions: translating the DNA sequence to an amino acid sequence, reading in the fasta files and reading in the blosum matrix. We tested these functions and added different blossum matrices. Additionally we added the main function and discussed the input paramteters and the output format. Next we started with the actual frameshift aware alginment. Therefore we decided to implement two algorithms in pairs of two and compared both implementations. This allowed us to have two different approaches where we could discuss difficulties and pair programming enabled efficient implementing of the algorithm. After the implementation of the frameshift aware alignment, we considered integrating the affine gap penalty. After some research we came up with a paper that seemed to solve this problem. However, they used a different approach doing the frameshift alignment than we did. Consequently, one group decided to implement the affine gap penalty by adding the affine gap penalty restrictions to our affine gap penalty implementation. The other group tried to implement as described in the paper. However, this approach seemed to be too challenging which is why we decided for the other implementation. Therefore working in pairs to resolve the tasks was a good idea.



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


