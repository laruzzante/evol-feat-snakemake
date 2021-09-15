# Compute orthogroup number of expansions across the tree (counts of gene number increases from the CAFE5 Base_change results table)

import pickle
from collections import defaultdict

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
copy_number_variation_table = snakemake.input.copy_number_variation_table
output_file_orthogroups = open(snakemake.output[0], 'w')

expansions = {}

with open(copy_number_variation_table) as f:
    for line in f:
        splitline = line.strip().split('\t')
        orthogroup = splitline[0]
        EXP = splitline[1]
        expansions[orthogroup] = EXP


# Process output files

output_file_orthogroups.write('orthogroup' + '\t' + 'EXP' + '\n')

for orthogroup in sorted(orthogroups.keys()):
    if orthogroup in expansions.keys():
        EXP = expansions[orthogroup]
        output_file_orthogroups.write(orthogroup + '\t' + str(EXP) + '\n')
    else:
        print(f"WARNING: {orthogroup} has no CAFE gene copy-number expansions (EXP) value.")


# Close files
output_file_orthogroups.close()
