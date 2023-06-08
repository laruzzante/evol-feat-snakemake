import sys

if len(sys.argv) > 1:
    infile = sys.argv[1]
else:
    print('Specify odb_OG2gene.tab input file.')


species_dict = {}

with open(infile, 'r') as f:
    for line in f:
        species_and_gene = line.strip().split('\t')[1]
        species_id = species_and_gene.split(':')[0]
        if species_id not in species_dict.keys():
            species_dict[species_id] = ''

print('\n'.join(species_dict.keys()))
