# Merge all computed features and print them in tsv file

import pandas as pd
from functools import reduce

# Retrieve information from Snakemake
orthogroup_features_files = snakemake.input
# input_files_genes = snakemake.input.gene_features
output_file = open(snakemake.output[0], 'w')
# output_file_genes = open(snakemake.output["gene_features"], 'w')

# Process output files
orthogroup_features = []
for input_file in orthogroup_features_files:
    df = pd.read_csv(open(input_file), sep='\t')
    # I <open(input_file)> and don't just use <input_file> here because if given a string,
    # pandas.read_csv() will automatically close the file, which is instead needed to be open for later.
    orthogroup_features.append(df)

merged_orthogroup_features = reduce(lambda left, right: pd.merge(left, right, on='orthogroup', how='outer'), orthogroup_features).fillna('NA')
merged_orthogroup_features.to_csv(output_file, sep='\t', index=False)

# gene_features = []
# for input_file in input_files_genes:
#     df = pd.read_csv(open(input_file), sep='\t')
#     gene_features.append(df)

# merged_gene_features = reduce(lambda left,right: pd.merge(left, right, on='gene', how='outer'), gene_features).fillna('NA')
# merged_gene_features.to_csv(output_file_genes, sep='\t', index=False)

# Close files
output_file.close()
# output_file_genes.close()
