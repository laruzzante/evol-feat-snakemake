#!/usr/bin/env python3

# Script that identifies orthogroups with high variation across species copy-number.
# Useful for filtering out problematic families in CAFE analysis.

import pickle
import statistics as stat

orthogroups = pickle.load(open('/home/lruzzant/evol-feat-snakemake/workflow/output/.orthogroups.pickle', 'rb'))

og2speccounts = {}

for og in orthogroups:
    for spec in set(orthogroups[og]["species"]):
        speccount = orthogroups[og]["species"].count(spec)
        if og not in og2speccounts.keys():
            og2speccounts[og] = [speccount]
        else:
            og2speccounts[og].append(speccount)

for og in sorted(og2speccounts):
    var = stat.variance(og2speccounts[og])
    print(f"{og}\t{var}")
