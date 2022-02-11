#!/usr/bin/env python3

'''
Author: Livio Ruzzante
Date: 27.01.2022
Script that takes an OrthoDB og2gene file with OrthoDB gene ids
and a GFF file with VectorBase gene ids, and creates a dicionary
map between the two gene id naming schemes while creating a flat
pickle python dictionary for it, for later usage in other python
scripts.
'''

import pickle

odbid_2_vbid = {}
odbid_2_vbid_mismatch = {}

vbid_2_odbid = {}
vbid_2_odbid_mismatch = {}

odbids = {}

with open('at6656_odb10v1_OG2genes.tab') as f:
    for line in f:
        odbid = line.strip().split('\t')[1]
        odbids[odbid] = None
    
print(len(odbids))

with open('odb10v1_genes.tab') as f2:
    for line in f2:
        odbid = line.strip().split('\t')[0]
        vbid = line.strip().split('\t')[2]
        if odbid in odbids.keys():
            if odbid not in odbid_2_vbid.keys():
                odbid_2_vbid[odbid] = vbid
            else:
                odbid_2_vbid_mismatch[odbid] = [odb_2_vbid[odbid]]
                odb_2_vbid_mismatch[odbid].append(vbid)
                print(odb_2_vbid_mismatch)
            if vbid not in vbid_2_odbid.keys():
                vbid_2_odbid[vbid] = odbid
            else:
                vbid_2_odbid_mismatch[vbid] = [vbid_2_odbid[vbid]]
                vbid_2_odbid_mismatch[vbid].append(odbid)
                print(vbid_2_odbid_mismatch)

print('odb_2_vbid_mismatches:')
print(len(odbid_2_vbid_mismatch))
print('vbid_2_odbid_mismtaches:')
print(len(vbid_2_odbid_mismatch))

print('Test prints')
print(len(odbid_2_vbid))
print(len(vbid_2_odbid))

print(odbids.popitem())
print(odbid_2_vbid.popitem())
print(vbid_2_odbid.popitem())

if len(odbid_2_vbid_mismatch) == 0 and len(vbid_2_odbid_mismatch) == 0:
    print('No mismatches between gene id naming schemes')
else:
    print('Mismatches present')


with open('vbid_2_odbid.pickle', 'wb') as handle:
    pickle.dump(vbid_2_odbid, handle, protocol=pickle.HIGHEST_PROTOCOL)

