# Create dictionaries of gene-orthogroup-species relationships from orthology table

import pandas as pd
import pickle

# Retrieve information from Snakemake
input_file = snakemake.input[0]
output_file_orthogroups = open(snakemake.output.orthogroups, 'wb') # 'wb' to write binary object instead of strings
output_file_genes = open(snakemake.output.genes, 'wb')
output_file_info = open(snakemake.output.info, 'w')

# I <open(input_file)> and don't just use <input_file> here because if given a string,
# pandas.read_csv() will automatically close the file, which is instead needed to be open for later.
df = pd.read_csv(input_file, sep='\t')
n_species = len(set(df['species']))
n_orthogroups = len(set(df['orthogroup']))
n_genes = len(set(df['gene']))

output_file_info.write('n_species\tn_orthogroups\tn_genes\n')
output_file_info.write(f'{n_species}\t{n_orthogroups}\t{n_genes}\n')

# Check consitency between number of species in config file and present in input orthology table
# n_species_from_config = float(snakemake.config["n_species"])
# if n_species_computed != n_species_from_config:
#     print('WARNING: number of species specified in config.yaml file different \
# from number of species found in <', snakemake.input[0] , '>.')

# Create dictionaries of orthogroups and gene features
orthogroups = {}
genes = {}

with open(input_file) as f:
    for line in f:
        line = line.strip()
        if 'orthogroup' in line or 'gene' in line or 'species' in line:
            continue
        orthogroup = line.split('\t')[0]
        gene = line.split('\t')[1]
        spec = line.split('\t')[2]
        if orthogroup not in orthogroups.keys():
            orthogroups[orthogroup] = {"genes": [gene], "species": [spec]}
        else:
            orthogroups[orthogroup]["genes"].append(gene)
            orthogroups[orthogroup]["species"].append(spec)
        if gene not in genes.keys():
            genes[gene] = {"orthogroup": [orthogroup], "species": [spec]}
        else:
            genes[gene]["orthogroup"].append(orthogroup)
            genes[gene]["species"].append(spec)

# Check integrity of orthology table
for gene in genes.keys():
    if len(genes[gene]["orthogroup"]) > 1:
        print('WARNING: gene', gene, 'mapped to more than 1 orthogroup.')
    if len(genes[gene]["species"]) > 1:
        print('WARNING: gene', gene, 'associated to more than 1 species.')

# Process output files
pickle.dump(orthogroups, output_file_orthogroups, protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(genes, output_file_genes, protocol=pickle.HIGHEST_PROTOCOL)

# Close files
output_file_orthogroups.close()
output_file_genes.close()
output_file_info.close()
