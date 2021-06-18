# Merge all computed features and print them in tsv file

import pandas as pd
from functools import reduce

# Retrieve information from Snakemake
gene_features_files = snakemake.input
output_file = open(snakemake.output[0], 'w')

# Process output files
gene_features = []
for input_file in gene_features_files:
    # I <open(input_file)> and don't just use <input_file> here because if given a string,
    # pandas.read_csv() will automatically close the file, which is instead needed to be open for later.
    df = pd.read_csv(open(input_file), sep='\t')
    gene_features.append(df)

if gene_features != []:
    merged_gene_features = reduce(lambda left, right: pd.merge(left, right, on='gene', how='outer'), gene_features).fillna('NA')
    merged_gene_features.to_csv(output_file, sep='\t', index=False)
else:
    print('No gene features.')

# Close files
output_file.close()
