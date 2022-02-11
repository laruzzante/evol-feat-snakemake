#!/usr/bin/env python3

'''
Author: Livio Ruzzante
Date: 27.01.2022
Script that takes a pickle dictionary of gene id mapping
and rewrites a GFF file replacing its gene ids with
OrthoDB gene ids.
'''

import sys
import pickle

gff = sys.argv[1]

print(gff)

with open('vbid_2_odbid.pickle', 'rb') as handle:
    vbid_2_odbid = pickle.load(handle)

converted_gff = gff + '_ODBid'

with open(gff) as f, open(converted_gff, 'w') as f2:
    for line in f:
        id = line.strip().split('\t')[8]
        start = id.find('=')
        end = id.find('-')
        id = id[start+1:end]
        newline = line.strip().split('\tID=')[0]
        newline = newline + '\t' + id + '\n'
        f2.write(newline)

