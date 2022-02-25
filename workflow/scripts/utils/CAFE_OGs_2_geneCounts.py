#!/usr/bin/env python3

import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Producing \
    a file with gene counts per OG per species from the OrthoDB \
    OGs default file (E.g. Insecta_OGs.txt)')

parser.add_argument('-i', '--input', dest='inputFile', metavar='input_file',
                    type=str, help='input file containg OGs at specific node.')

args = parser.parse_args()

outputName = args.inputFile.replace('.tab', '_geneCounts.tab')

with open(args.inputFile) as f:
    lines = f.readlines()

ogDict = defaultdict(lambda: defaultdict(list))
ogCountsDict = defaultdict(dict)
specList = []

for line in lines[1:]:

    line = line.strip()
    og = line.split('\t')[0]
    geneid = line.split('\t')[1]
    spec = line.split('\t')[2]

    ogDict[og][spec].append(geneid)
    specList.append(spec)

all_species = sorted(set(specList))
all_species_str = '\t'.join(all_species)

print('OG_species_span\tOG\t' + all_species_str)

for og in sorted(ogDict.keys()):
    for spec in all_species:
        ogCountsDict[og][spec] = 0

    for specInOG, genes in sorted(ogDict[og].items()):
        ogCountsDict[og][specInOG] += len(genes)

for og, species in sorted(ogCountsDict.items()):
    outline = str(len(ogDict[og].keys())) + '\t' + og
    for spec in sorted(species):
        outline = outline + '\t' + str(ogCountsDict[og][spec])
    print(outline)
