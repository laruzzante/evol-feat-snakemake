#!/usr/bin/env python3
import sys
from collections import defaultdict
from shutil import copyfile
import os
import subprocess

## This script removes alternative transcripts protein sequences in a gene sets file. Keeps the longest protein sequence.
## Useful when assessing BUSCOs in -prot mode.

## Check if the required input file is passed to the script
if( len(sys.argv) > 1 ):
    infile = open(sys.argv[1])
else:
    print('Please specify gene set input file e.g. python3 filterGeneSets_BUSCO.py input_FASTA_file')
    sys.exit()

## Counting gene duplicates (differentially expressed) and appending sequences to dict.
protDic = defaultdict(list)
dupliDic = defaultdict(list)
seqDic = defaultdict(list)
isSeq = False
lines = infile.readlines()
l = 0

for line in lines:
    l += 1

    if( line.startswith('>') ):
        if isSeq:
            seqDic[prot].append(sequence)
        isSeq = False
        sequence = str()
        prot = line.split('\n')[0].strip('\n')
        prot = prot.strip('>')
        prot = prot.split(' ')[0]

        ## If the protein ID has the form AB_123456.1, it only keeps the integer part of the number,
        ## as the decimals are possibly alternative transcripts of the same gene.
        if( '.' in prot ):
            prot = prot.split('.')[0]

        if( prot not in protDic.keys() ):
            protDic[prot] = 1

        else:
            dupliDic[prot] = ''
            protDic[prot] += 1
    else:
        sequence = sequence + line.strip('\n')
        isSeq = True

    if( l == len(lines) ):
        seqDic[prot].append(sequence)


print('\n')

for k,v in protDic.items():
    if( protDic[k] > 1 ):
        print(k+' is present '+str(v)+' times in the given gene set.')

if( not bool(dupliDic) ):
    print('No differentially expressed genes found in the given gene set.')

## Writing sequences to new file, without duplicates.
outfileName = sys.argv[1].split('.faa')[0]
outfileName = outfileName + '_noAltTranscripts.faa'

print('\nFile with only the longest protein sequences from differentially expressed genes saved as < '+outfileName+' >.\n')

with( open(outfileName,'w') ) as outfile:
    for k,v in seqDic.items():
        if( k not in dupliDic ):
            outfile.write('>'+str(k)+'\n')
            outfile.write(str(v).strip('[').strip(']').strip('\'')+'\n')
        else:
            outfile.write('>'+str(k)+'\n')
            ## Writing only longest protein sequence of the differentially expressed genes.
            outfile.write(max(seqDic[k],key=len)+'\n')
