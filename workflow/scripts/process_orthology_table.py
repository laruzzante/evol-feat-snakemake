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

    # Creating dictionaries from formatted orthology table
    for line in f:
        line = line.strip()
        if 'orthogroup' in line or 'gene' in line or 'species' in line:
            continue

        # Assuming the required format is satisfied
        orthogroup = line.split('\t')[0]
        gene = line.split('\t')[1]
        if len(line.split('\t')) > 2:
            spec = line.split('\t')[2]
        # However if column species is not provided, we try to infer gene
        # and species name from OrthoDB gene codes.
        elif ':' in gene:
            spec = gene.split(':')[0]
        else:
            print('WARNING: no species column detected. Please format your orthology table as required.')

        # Adding orthogroup to dictionary
        if orthogroup not in orthogroups.keys():
            # Adding the species in this way will provide a list of species to each orthogroup
            # which can be repeated. I.E. multiple appearances of the same species name
            # might appear in the orthogroup species list, due to multicopy genes, hence to
            # get a unique species list we will still have to "set" the list. Leaving it like
            # this is useful to compute later multi-copy evolutionary features.
        else:
            orthogroups[orthogroup]["genes"].append(gene)
            orthogroups[orthogroup]["species"].append(spec)
        if gene not in genes.keys():
            genes[gene] = {"orthogroup": orthogroup, "species": spec}
        else:
            with open(snakemake.log.log, "a") as logfile:
                logfile.write(f'WARNING: gene {gene} already associated to orthogroup \
{genes[gene]["orthogroup"]}. Association with {orthogroup} will be discarded.\n')
            # genes[gene]["orthogroup"].append(orthogroup)
            # genes[gene]["species"].append(spec)
        if spec not in species.keys():
            species[spec] = {"orthogroups": [orthogroup], "genes": [gene]}
        else:
            species[spec]["orthogroups"].append(orthogroup)
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
# Saving dictionaries as pickle files so that they can be easily used by other python scripts later on.
pickle.dump(orthogroups, output_file_orthogroups, protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(genes, output_file_genes, protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(species, output_file_species, protocol=pickle.HIGHEST_PROTOCOL)

# Close files
output_file_orthogroups.close()
output_file_genes.close()
output_file_info.close()
