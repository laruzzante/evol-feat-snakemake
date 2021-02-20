import os
import pandas as pd

filepaths = snakemake.input
output_dir = snakemake.output[0]
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

for input_file in filepaths:
    df = pd.read_csv(input_file, sep='\t')
    gene_column_header = df.columns[0]
    for feature in df.columns[1:]:
        df_subset = df[[gene_column_header, feature]]
        df_subset = df.rename(columns={gene_column_header: 'gene'})
        output_file = os.path.join(output_dir, f'{feature}.tsv')
        df_subset.to_csv(output_file, sep='\t', index=False, na_rep='NA')
