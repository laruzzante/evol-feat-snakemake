# Merge all computed features and print them in tsv file

import pandas as pd
from functools import reduce

# Retrieve information from Snakemake
input_files_orthogroups = snakemake.input["evol_feat_orthogroups"]
input_files_genes = snakemake.input["evol_feat_genes"]
output_file_orthogroups = open(snakemake.output["evol_feat_orthogroups"], 'w')
output_file_genes = open(snakemake.output["evol_feat_genes"], 'w')

# Process output files
orthogroup_features = []
for input_file in input_files_orthogroups:
    df = pd.read_csv(open(input_file), sep='\t') # I don't use the input_file variable name here because if given a string, pandas.read_csv() will automatically close th file, which is instead needed to be open for later
    orthogroup_features.append(df)

merged_orthogroup_features = reduce(lambda left, right: pd.merge(left, right, on='orthogroup', how='outer'), orthogroup_features).fillna('NA')
merged_orthogroup_features.to_csv(output_file_orthogroups, sep='\t', index=False)

gene_features = []
for input_file in input_files_genes:
    df = pd.read_csv(open(input_file), sep='\t')
    gene_features.append(df)

merged_gene_features = reduce(lambda left,right: pd.merge(left, right, on='gene', how='outer'), gene_features).fillna('NA')
merged_gene_features.to_csv(output_file_genes, sep='\t', index=False)

# Close files
output_file_orthogroups.close()
output_file_genes.close()
