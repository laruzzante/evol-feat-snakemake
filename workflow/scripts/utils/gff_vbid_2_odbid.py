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
gff_converted = sys.argv[2]

print(gff)

with open('vbid_2_odbid.pickle', 'rb') as handle:
    vbid_2_odbid = pickle.load(handle)


with open(gff) as f, open(gff_converted, 'w') as f2:
    for line in f:
        id = line.strip().split('\t')[8]
        source = line.strip().split('\t')[1]

        case0 = ['ID=', '-']
        case1 = ['ID=cds-',';']
        case2 = ['Name=', '-']
        case3 = ['Name=', ';']
        case4 = ['Parent=', '-']
        case5 = ['Parent=', '.t1']
        case6 = ['Parent=', '.t2']

        if id.find(case0[0]) != -1 and id.find(case0[1]) != -1 and id.find('=cds') == -1 and id.find('=CDS') == -1:
            start = id.find(case0[0])
            end = id.find(case0[1])
            case = case0
        elif id.find(case1[0]) != -1 and id.find(case1[1]) != -1:
            start = id.find(case1[0])
            end = id.find(case1[0])
            case = case1
        elif id.find(case2[0]) != -1 and id.find(case2[1]) != -1:
            start = id.find(case2[0])
            end = id.find(case2[0])
            case = case2
        elif id.find(case3[0]) != -1 and id.find(case3[1]) != -1:
            start = id.find(case3[0])
            end = id.find(case3[0])
            case = case3
        elif id.find(case4[0]) != -1 and id.find(case4[1]) != -1:
            start = id.find(case4[0])
            end = id.find(case4[0])
            case = case4
        elif id.find(case5[0]) != -1 and id.find(case5[1]) != -1:
            start = id.find(case5[0])
            end = id.find(case5[0])
            case = case5
        elif id.find(case6[0]) != -1 and id.find(case6[1]) != -1:
            start = id.find(case6[0])
            end = id.find(case6[0])
            case = case6
        else:
            print('WARNING: could not extract gene name in:')
            print(id)

        # Specific cases that I could not solve with a generic rule
        if source == 'Genbank':
            case = ['ID=cds-', ';']
            start = id.find(case[0])
            end = id.find(case[1])
        elif source == 'OGS1':
            case = ['Name=', '-']
            start = id.find(case[0])
            end = id.find(case[1])
        elif source == 'FlyBase':
            case = ['Parent=', '']
            start = id.find(case[0])
            end = len(id) + 1
        elif source == 'EMBL':
            case = ['ID=cds-', ';']
            start = id.find(case[0])
            end = id.find(case[1])
        elif source == 'RefSeq' or source == 'AUGUSTUS' \
            or source == 'bimp_OGSv1.0' or source == 'Melipona_quadrifasciata_v1.1_HGD' \
            or source == 'Labl_OGS_v5.42':
            case = ['Parent=', '-']
            start = id.find(case[0])
            if id.find('.t') != -1 and id.find('-') == -1:
                case = ['Parent=', '']
                end = len(id) + 1
            else:
                end = id.find(case[1])

        try:
            id = id[start+len(case[0]):end+len(case[1])-1]
            newline = line.strip().split(case[0])[0]
            newline = newline.rsplit('\t', 1)[0]
            newline = newline + '\t' + id + '\n'
            f2.write(newline)
        except NameError:
            continue
