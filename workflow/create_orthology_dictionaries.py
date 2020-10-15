# Compute orthogroup universality (relative species-span)

import pandas as pd
import pickle

# Retrieve information from Snakemake
input_file = open(snakemake.input[0])
output_file_orthogroups = open(snakemake.output["orthogroups"], 'wb')
output_file_genes = open(snakemake.output["genes"], 'wb')

# Check consitency between number of species in config file and present in input orthology table
df = pd.read_csv(open(snakemake.input[0]), sep='\t') # I don't use the input_file variable name here because if given a string, pandas.read_csv() will automatically close th file, which is instead needed to be open for later
n_species_computed = len(set(df['species']))
n_species_from_config = float(snakemake.config["n_species"])
if n_species_computed != n_species_from_config:
    print('WARNING: number of species specified in config.yaml file different \
from number of species found in <', snakemake.input[0] , '>.')

# Create dictionaries of orthogroups and gene features
orthogroups = {}
genes = {}
for line in input_file:
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
input_file.close()
output_file_orthogroups.close()
output_file_genes.close()
