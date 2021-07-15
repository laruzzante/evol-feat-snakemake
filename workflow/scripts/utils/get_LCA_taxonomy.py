#!/usr/bin/env python3

from collections import defaultdict
from ete3 import NCBITaxa
import sys

ncbi = NCBITaxa()
ncbi.update_taxonomy_database()


if len(sys.argv) > 1:
    infile = open(sys.argv[1], 'r')
else:
    print('Specify FASTA input file.')

odb = infile.readlines()
infile.close()

ogdict = defaultdict(list)

print('og_id' + '\t' + 'last_common_ancestor' + '\t' + 'species_list')

for line in odb:
    og = line.strip().split('\t')[0]
    gene = line.strip().split('\t')[1]
    taxid = gene.split('_')[0]

# Managing exceptions, where NCBI taxids were updated after the OrthoDB10 update
    if taxid == '1415176':
        taxid = '2587831'

    ogdict[og].append(taxid)

for og in ogdict.keys():
    lineages = []
    commonAncestors = []
    isCommonAncestor = True
    for taxid in set(ogdict[og]):
        lineage = ncbi.get_lineage(taxid)
        lineages.append(lineage)

    for rank in lineages[0]:
        for lineage in lineages[1:]:
            if rank not in lineage:
                isCommonAncestor = False
                break
        if isCommonAncestor is True:
            commonAncestors.append(rank)

    try:
        commonAncestors_names = ncbi.get_taxid_translator(commonAncestors)
        LCA = commonAncestors[-1]
        LCA_name = commonAncestors_names[LCA]

        species = ncbi.get_taxid_translator(set(ogdict[og]))
        species = sorted(list(species.values()))

        print(og + '\t' + LCA_name + '\t' + ', '.join(species) + '\n')

    except ValueError:
        break
