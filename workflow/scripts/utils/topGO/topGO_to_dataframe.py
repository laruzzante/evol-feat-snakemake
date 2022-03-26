#!/usr/bin/env python3
from collections import defaultdict

# Parsing AGAMP GO terms
with open("AGAP2GOS-topGO.txt") as f:
    GOlines = f.readlines()

GOdict = defaultdict(list)
for line in GOlines:
    splitline = line.strip().split('\t')
    geneid = splitline[0]
    GOterms = splitline[1].split(', ')
    for GOterm in GOterms:
        GOdict[geneid].append(GOterm)
        print(geneid + '\t' + GOterm)
