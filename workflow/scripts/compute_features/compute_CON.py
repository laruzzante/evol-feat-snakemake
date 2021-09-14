# Compute orthogroup number of contractions across the tree (counts of gene number losses from the CAFE5 Base_change results table)

import pickle
from collections import defaultdict

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
copy_number_variation_table = snakemake.input.copy_number_variation_table
output_file_orthogroups = open(snakemake.output[0], 'w')

contractions = {}

with open(copy_number_variation_table) as f:
    for line in f:
        splitline = line.strip().split('\t')
        orthogroup = splitline[0]
        CON = splitline[3]
        contractions[orthogroup] = CON


# Process output files

output_file_orthogroups.write('orthogroup' + '\t' 'CON' + '\n')

for orthogroup in sorted(orthogroups.keys()):
    if orthogroup in contractions.keys():
        CON = contractions[orthogroup]
        output_file_orthogroups.write(orthogroup + '\t' + str(CON) + '\n')
    else:
        print(f"WARNING: {orthogroup} has no CAFE gene copy-number contractions (CON) value.")


# Close files
output_file_orthogroups.close()
