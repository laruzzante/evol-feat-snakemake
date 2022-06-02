import sys
import os

if len(sys.argv) > 1:
    infile = open(sys.argv[1], 'r')
else:
    print('Specify input file.')
    sys.exit()

outdir = 'GOFigure_inputs_' + sys.argv[1]
outdir = outdir.strip('.tsv')
os.mkdir(outdir)

for line in infile:
    if line.startswith('Cluster'):
        cluster = line.strip().split(' ')[1]
        filename = cluster + '.tsv'
        outpath = os.path.join(outdir, filename)
        continue
    elif line.startswith('GO.ID'):
        continue
    else:
        goterm = line.split('\t')[0]
        pval = line.split('\t')[7].strip()
        with open(outpath, 'a') as outfile:
            outfile.write(goterm + '\t' + pval + '\n')

infile.close()
