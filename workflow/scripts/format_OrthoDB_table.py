import pandas as pd

orthology_table_path = snakemake.input[0]
output_file = open(snakemake.output[0], 'w')

# I <open(input_file)> and don't just use <input_file> here because if given a string,
# pandas.read_csv() will automatically close the file, which is instead needed to be open for later.
df = pd.read_csv(open(orthology_table_path), sep='\t')

needs_columns_formatting = False
needs_geneid_formatting = False

# Checking that the first column has the format seen from OrthoDB output
if df.columns[0] == 'pub_og_id':
    needs_columns_formatting = True
    # Checking that gene identifiers are unqiue
    if not df['pub_gene_id'].is_unique:
        needs_geneid_formatting = True

# If gene identifiers are not unique, add species code in front of the gene id
if needs_geneid_formatting:
    df['pub_gene_id'] = df['code'].astype(str) + ':' + df['pub_gene_id'].astype(str)

# Extracting only relevant columns and rename
if needs_columns_formatting:
    df = pd.DataFrame(df[['pub_og_id', 'pub_gene_id', 'code']])
    df = df.rename(columns={'pub_og_id': 'orthogroup', 'pub_gene_id': 'gene', 'code': 'species'})

# Saving as formatted orthoology table
df.to_csv(output_file, sep='\t', index=False)

# Close files
output_file.close()
