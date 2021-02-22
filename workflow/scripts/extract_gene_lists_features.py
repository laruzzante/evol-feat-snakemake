import  pandas as pd

gene_list_files = snakemake.input
merged_genes_file = snakemake.input[1]
output_file = open(snakemake.output[0], 'w')

merged_genes_df = pd.read_csv(merged_genes_file, sep='\t')

for file in gene_list_files:
    gene_names = []
    with open(file) as f:
        for line in f:
            gene_names.append(line.strip())
        # gene_names = f.readlines()
    # gene_names.strip()
    # print(merged_genes_df['gene'][0:10])
    # print(gene_names)
    df = merged_genes_df.loc[merged_genes_df['gene'].isin(gene_names)]
        # for line in f[:-1]:
            # if line in merged_genes_df['gene']:
                # outfile.write(merged_genes_df[])
    # for
df.to_csv(output_file, sep='\t', index=False)

# Close files
output_file.close()
