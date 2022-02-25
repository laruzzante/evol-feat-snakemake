#!/usr/bin/env python3

'''
Author: Livio Ruzzante
Date: 22.02.2022
Script that takes a pickle dictionary of gene id mapping
and rewrites a GFF file replacing its gene ids with
OrthoDB gene ids.
'''

import os
import pickle

vbid_2_odbid = pickle.load(open('/media/lruzzant/DATA/evol-feat-data-processing/vbid_2_odbid.pickle', 'rb'))
odbid_2_vbid = pickle.load(open('/media/lruzzant/DATA/evol-feat-data-processing/odbid_2_vbid.pickle', 'rb'))


directory = os.fsencode('/media/lruzzant/DATA/evol-feat-data-processing/GFFs_GENES/')
out_directory = os.fsencode('/media/lruzzant/DATA/evol-feat-data-processing/GFFs_GENES_ODB/')

for file in os.listdir(directory):
    gff_file = os.fsdecode(file)
    outfile = gff_file + '_ODB'
    taxid = gff_file.split('_')[0] + '_0'
    print('Processing :' + outfile)
    print('\t- taxid: ' + taxid)

    odb_gene_ids = {}
    for key in odbid_2_vbid.keys():
        if taxid in key:
            odb_gene_ids[key] = odbid_2_vbid[key]

    num_lines = sum(1 for line in open(os.path.join(os.fsdecode(directory), gff_file)))
    i = 0

    with open(os.path.join(os.fsdecode(directory), gff_file)) as gff,\
     open(os.path.join(os.fsdecode(out_directory), outfile), 'w') as outfile:
        for line in gff:
            i += 1
            if i%1000 == 0:
                progress = round(100*i/num_lines , 2)
                print('\t\t' + str(progress) + ' %')
            col8 = line.strip().split('\t')[8]
            for key in odb_gene_ids.keys():
                if odb_gene_ids[key] in col8:
                    newline = line.rsplit('\t', 1)[0]
                    newline = newline + '\t' + key + '\n'
                    outfile.write(newline)
