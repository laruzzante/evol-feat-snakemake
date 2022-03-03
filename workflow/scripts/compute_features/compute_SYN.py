# Compute orthogroup synteny

import pickle

# Retrieve information from Snakemake
genes = pickle.load(open(snakemake.input.genes, 'rb'))
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
ordered_gff_genes = pickle.load(open(snakemake.input.ordered_gff_genes, 'rb'))
output_file_orthogroups = open(snakemake.output[0], 'w')
