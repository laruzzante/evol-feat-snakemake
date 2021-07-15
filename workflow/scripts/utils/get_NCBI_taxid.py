#!/usr/bin/env python3
import sys
from ete3 import NCBITaxa

ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()

print('\nTAXID usage example: ./taxid.py Tribolium castaneum')
print('TAXID output example: 7070, Tribolium castaneum;\n')

taxnames = []
taxname = ' '.join(sys.argv[1:])
taxnames.append(taxname)

try:
    taxids = ncbi.get_name_translator(taxnames)
    taxid = str(taxids[taxname][0])
except TypeError:
    print('Input TypeError. Input only takes a species name.\nUsage example: ./taxid.py Apis mellifera')

print(taxid + ',', taxname)
