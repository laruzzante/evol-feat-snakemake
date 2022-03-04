# Compute orthogroup number of stabilities across the tree (counts of gene number maintenance from the CAFE5 Base_change results table)

import pickle

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
copy_number_variation_table = snakemake.input.copy_number_variation_table
output_file_orthogroups = open(snakemake.output[0], 'w')

stabilities = {}

with open(copy_number_variation_table) as f:
    next(f)
    for line in f:
        splitline = line.strip().split('\t')
        orthogroup = splitline[0]
        STA = splitline[2]
        stabilities[orthogroup] = STA


# Process output files

output_file_orthogroups.write('orthogroup' + '\t' + 'STA' + '\n')

for orthogroup in sorted(orthogroups.keys()):
    if orthogroup in stabilities.keys():
        STA = stabilities[orthogroup]
        output_file_orthogroups.write(orthogroup + '\t' + str(STA) + '\n')
    else:
        print(f"WARNING: {orthogroup} has no CAFE gene copy-number stabilities (STA) value.")


# Close files
output_file_orthogroups.close()
