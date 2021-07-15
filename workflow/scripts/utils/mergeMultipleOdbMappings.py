#!/usr/bin/env python3

import sys
from collections import defaultdict


# Script that merges 2 geneset mapping to ODB files.

if (len(sys.argv) == 3):
    map1 = open(sys.argv[1], 'r')
    map2 = open(sys.argv[2], 'r')
else:
    print('Error: specify 2 ODB mapping files.')
    print('Usage: ./mergeMultipleOdbMappings.py MappingFile1 MappingFile2')
    sys.exit()

map1lines = map1.readlines()
map2lines = map2.readlines()

headerLine = map1lines[0]
map1LinesDic = defaultdict(dict)
map2LinesDic = defaultdict(dict)
mapMergedLinesDics = defaultdict(dict)

for line in map1lines[1:]:
    ogID = line.strip('\n').split(' ')[0]
    geneID = line.strip('\n').split(' ')[1]
    map1LinesDic[ogID][geneID] = line

for line in map2lines[1:]:
    ogID = line.strip('\n').split(' ')[0]
    geneID = line.strip('\n').split(' ')[1]
    map2LinesDic[ogID][geneID] = line

mapMergedLinesDics = map1LinesDic


for ogID in sorted(map2LinesDic.keys()):
    for geneID in sorted(map2LinesDic[ogID].keys()):
        if (ogID not in sorted(map1LinesDic.keys())):
            mapMergedLinesDics[ogID][geneID] = map2LinesDic[ogID][geneID]
        elif (ogID in sorted(map1LinesDic.keys()) and geneID not in sorted(map1LinesDic[ogID].keys())):
            mapMergedLinesDics[ogID][geneID] = map2LinesDic[ogID][geneID]


# outfileName = "mergedMaps.txt"
# outfile = open(outfileName, 'w')
# outfile.write(headerLine)
print(headerLine.strip())
for ogID in sorted(mapMergedLinesDics.keys()):
    for geneID in sorted(mapMergedLinesDics[ogID].keys()):
        print(mapMergedLinesDics[ogID][geneID].strip())
        # outfile.write(mapMergedLinesDics[ogID][geneID])
# outfile.close()
# print("File saved as", outfileName)
