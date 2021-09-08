import pandas as pd
import sys

orthology_table_path = snakemake.input.orthology_table
output_file = open(snakemake.output.formatted_orthology_table, 'w')

isOrthoDBv9 = False
isOrthoDBv10 = False
hasHeader = False
hasSpeciesCol = False

# Checking if the provided orthology table has the required default header
default_header = 'orthogroup\tgene\tspecies'

with open(orthology_table_path, 'r') as infile:
    try:
        first_line = infile.readline()
    except EOFError:
        print('EOFerror: empty orthology input file.')

    if default_header in first_line:
        hasHeader = True
    else:
        print('WARNING: cannot retrieve required header from the orthology table.')
        print('\tThe suggested header format is: orthogroup <tab> gene <tab> species')
        print('\tAssuming no header is present, continuing to usual formatting...')

    if '\t' in first_line:
        splitted_first_line = first_line.split('\t')
        if len(splitted_first_line) > 2:
            hasSpeciesCol = True
            if 'pub_og_id' in splitted_first_line[0]:
                isOrthoDBv9 = True
            if len(splitted_first_line) > 3:
                print('WARNING: first row (header or data) has more than 3 columns.')
                print('\tColumns after Column #3 will be ignored.')
        elif len(splitted_first_line) == 2:
            if hasHeader == False and ':' in splitted_first_line[1]:
                isOrthoDBv10 = True
        else:
            print('ERROR: no tab-delimited data (or header), please format the orthology table as specified in the config file.')
            print('\tTerminating Evol-Feat.')
            sys.exit()


# I <open(input_file)> and don't just use <input_file> here because if given a string,
# pandas.read_csv() will automatically close the file, which is instead needed to be open for later.
def reading_OrthoDBv9():
    df = pd.read_csv(open(orthology_table_path), sep='\t')

    needs_columns_formatting = False
    needs_geneid_formatting = False

    ## ORTHODBv9
    # Checking that the first column has the format seen from OrthoDBv9 output
    if df.columns[0] == 'pub_og_id':
        needs_columns_formatting = True
    # Checking that gene identifiers are unqiue
    if df.columns[1] == 'pub_gene_id':
        if not df['pub_gene_id'].is_unique:
            needs_geneid_formatting = True
    else:
        'WARNING: gene id column not recognizable as OrthoDBv9 formatting'

    # If gene identifiers are not unique, add species code in front of the gene id
    if needs_geneid_formatting:
        df['pub_gene_id'] = df['code'].astype(str) + ':' + df['pub_gene_id'].astype(str)

    # Extracting only relevant columns and rename
    if needs_columns_formatting:
        df = pd.DataFrame(df[['pub_og_id', 'pub_gene_id', 'code']])
        df = df.rename(columns={'pub_og_id': 'orthogroup', 'pub_gene_id': 'gene', 'code': 'species'})

    return df

## ORTHODBv10
# In OrthoDB10 the orthology table has no headers and onyl 2 columns: orthogroup, geneid
# The geneid contains the species code (NCBI code followed by _X, usually 0 if arthropods)
# The _0 after the species code was implemented in OrthoDB in order to distinguish between
# different subspecies or viral strains, etc..
# Hence the geneid format looks like: NCBIspeciescode_0:geneid
def reading_OrthoDBv10():
    df = pd.read_csv(open(orthology_table_path), sep='\t', header=None)
    df.columns = ['orthogroup', 'gene']
    split_gene_names = df.gene.str.split(":", expand=True,)
    df['species'] = split_gene_names[0]

    return df


def wrong_format_message():
    print('ERROR: cannot retrieve species information.')
    print('\tSpecies column must be provided in Column #3 in tab-delimited rows.')
    print('\tSpecies names can also be automatically retrieved if present before <:> in the gene id column, as in default OrthoDBv10 format.')
    print('\tPlease reformat the input orthology table as specified in the config file.')
    print('\tTerminating Evol-Feat.')


if isOrthoDBv9:
    df = reading_OrthoDBv9()
    print('Assuming OrthoDBv9 input orthology table format.')
elif isOrthoDBv10:
    df = reading_OrthoDBv10()
    print('Assuming OrthoDBv10 input orthology table format.')
else:
    if hasHeader:
        if hasSpeciesCol:
            df = pd.read_csv(open(orthology_table_path), sep='\t')
        else:
            wrong_format_message()
            sys.exit()
    else:
        if hasSpeciesCol:
            df = reading_OrthoDBv10()
        else:
            df = wrong_format_message()
            sys.exit()

# Saving as formatted orthoology table
df.to_csv(output_file, sep='\t', index=False)

# Close files
output_file.close()
