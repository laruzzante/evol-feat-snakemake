# Computing orthogroup universality (relative species-span)

from collections import defaultdict
import pandas as pd

# Retrieve information from Snakemake
input_file = open(snakemake.input[0])
output_file_orthogroups = open(snakemake.output["UNI_orthogroups"], 'w')
output_file_genes = open(snakemake.output["UNI_genes"], 'w')


# Create dictionaries of orthogroups and gene features
orthogroups = defaultdict()
genes = defaultdict()

# df = pd.read_csv(input_file, sep='\t')
n_species = float(snakemake.config["n_species"]) #len(set(df['species']))

for line in input_file:
    orthogroup = line.split('\t')[0]
    spec = line.split('\t')[2]
    orthogroups[orthogroup]['species'].append(spec)
#
# # Process files
# for orthogroup in sorted(orthogroups.keys()):
#     output_file_orthogroups.write(orthogroup + orthogroup['species'] + '\n')
