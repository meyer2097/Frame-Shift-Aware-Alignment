# Frame-Shift-Aware-Alignment
Sequence Bioinformatics Group Project WS 21/22
Samira Breitling, Jules Kreuer, Sebastian Meyer and Catalina Schlotterer


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


