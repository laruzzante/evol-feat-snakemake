# Create dictionaries of gene-orthogroup-species relationships from orthology table

import pandas as pd
import pickle

# Retrieve information from Snakemake
input_file = snakemake.input[0]
output_file_orthogroups = open(snakemake.output.orthogroups, 'wb') # 'wb' to write binary object instead of strings
output_file_genes = open(snakemake.output.genes, 'wb')
output_file_species = open(snakemake.output.species, 'wb')
output_file_info = open(snakemake.output.info, 'w')


# Create dictionaries of orthogroups and gene features
orthogroups = {}
genes = {}
species = {}

with open(input_file) as f:
    for line in f:
        line = line.strip()
        if 'orthogroup' in line or 'gene' in line or 'species' in line:
            continue
        orthogroup = line.split('\t')[0]
        gene = line.split('\t')[1]
        spec = line.split('\t')[2]
        gene = spec + ':' + gene
        if orthogroup not in orthogroups.keys():
            orthogroups[orthogroup] = {"genes": [gene], "species": [spec]}
        else:
            orthogroups[orthogroup]["genes"].append(gene)
            orthogroups[orthogroup]["species"].append(spec)
        if gene not in genes.keys():
            genes[gene] = {"orthogroup": orthogroup, "species": spec}
        else:
            print(f'WARNING: gene {gene} already associated to orthogroup \
{genes[gene]["orthogroup"]}. Association with {orthogroup} will be discarded.')
            # genes[gene]["orthogroup"].append(orthogroup)
            # genes[gene]["species"].append(spec)
        if spec not in species.keys():
            species[spec] = {"orthogroup": [orthogroup], "genes": [gene]}
        else:
            species[spec]["orthogroup"].append(orthogroup)
            species[spec]["genes"].append(gene)

# Check integrity of orthology table
# for gene in genes.keys():
#     # Checking for association of gene id to more than 1 orthogroup
#     if len(genes[gene]["orthogroup"]) > 1:
#         # Check for multiple association of gene to same orthogroup
#         if len(genes[gene]["orthogroup"]) != len(set(genes[gene]["orthogroup"])):
#             # Removing duplicated orthogroup id
#             genes[gene]["orthogroup"] = set(genes[gene]["orthogroup"])
#         # If after removing duplicates, the gene is still associated to multiple orthogroups, print Warning
#         if len(genes[gene]["orthogroup"]) > 1:
#             print('WARNING: gene id', gene, 'associated to', len(genes[gene]["orthogroup"]),'orthogroups:', genes[gene]["orthogroup"])
#     # Checking for association of gene id to more than 1 species
#     if len(genes[gene]["species"]) > 1:
#         # Check for multiple association of gene to same species
#         if len(genes[gene]["species"]) != len(set(genes[gene]["species"])):
#             # Removing duplicated species id
#             genes[gene]["species"] = set(genes[gene]["species"])
#         # If after removing duplicates, the gene is still associated to multiple species, print Warning
#         if len(genes[gene]["species"]) > 1:
#             print('WARNING: gene id', gene, 'associated to', len(genes[gene]["species"]), 'species:', genes[gene]["species"])

# Writing out orthology table Info statistics file
n_species = len(species.keys())
n_orthogroups = len(orthogroups.keys())
n_genes = len(genes.keys())

output_file_info.write('n_species\tn_orthogroups\tn_genes\n')
output_file_info.write(f'{n_species}\t{n_orthogroups}\t{n_genes}\n')

# Process output files
pickle.dump(orthogroups, output_file_orthogroups, protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(genes, output_file_genes, protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(species, output_file_species, protocol=pickle.HIGHEST_PROTOCOL)

# Close files
output_file_orthogroups.close()
output_file_genes.close()
output_file_info.close()
