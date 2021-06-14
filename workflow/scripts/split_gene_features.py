import os
import pandas as pd

filepaths = snakemake.input
output_dir = snakemake.output[0]
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

for input_file in filepaths:
    df = pd.read_csv(input_file, sep='\t')
    gene_column_header = df.columns[0]

    # Because later on we have to merge these features files using the shared gene name column,
    # we must be sure that no gene name is present more than once in each feature file.
    # The gene names with multiple instances will be discarded from analysis with a warning.

    repeated_genes = df.duplicated(subset=gene_column_header, keep=False)
    if True in repeated_genes:
        print(f'WARNING: the {input_file} list contains repeated gene names.\n\n\
The following genes will be removed due to incompabilities with later feature merging rules:\n')
        print(', '.join(set(df[repeated_genes][gene_column_header])))
        df = df.drop_duplicates(subset=gene_column_header, keep=False, ignore_index=True)

    for feature in df.columns[1:]:
        df_subset = df[[gene_column_header, feature]]
        df_subset = df_subset.rename(columns={gene_column_header: 'gene'})
        output_file = os.path.join(output_dir, f'{feature}.tsv')
        df_subset.to_csv(output_file, sep='\t', index=False, na_rep='NA')
