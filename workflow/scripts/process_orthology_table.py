# Create dictionaries of gene-orthogroup-species relationships from orthology table

import pandas as pd
import pickle

# Retrieve information from Snakemake
input_file = snakemake.input.formatted_orthology_table
output_file_orthogroups = open(snakemake.output.orthogroups, 'wb') # 'wb' to write binary object instead of strings
output_file_genes = open(snakemake.output.genes, 'wb')
output_file_species = open(snakemake.output.species, 'wb')
output_file_orthogroups_2_species_2_genes = open(snakemake.output.orthogroups_2_species_2_genes, 'wb')
output_file_info = open(snakemake.output.info, 'w')


# Create dictionaries of orthogroups and gene features
orthogroups = {}
genes = {}
species = {}
orthogroups_2_species_2_genes = {}

with open(input_file) as f:
    next(f) # skipping header row of formatted table

    # Creating dictionaries from formatted orthology table
    for line in f:
        line = line.strip()

        # Assuming the required format is satisfied
        orthogroup = line.split('\t')[0]
        gene = line.split('\t')[1]
        spec = line.split('\t')[2]

        # Adding orthogroup to dictionaries
        if orthogroup not in orthogroups.keys():
            # Adding the species in this way will provide a list of species to each orthogroup
            # which can be repeated. I.E. multiple appearances of the same species name
            # might appear in the orthogroup species list, due to multicopy genes, hence to
            # get a unique species list we will still have to "set" the list. Leaving it like
            # this is useful to compute later multi-copy evolutionary features.
            orthogroups[orthogroup] = {"genes": [gene], "species": [spec]}
        else:
            orthogroups[orthogroup]["genes"].append(gene)
            orthogroups[orthogroup]["species"].append(spec)

        if gene not in genes.keys():
            genes[gene] = {"orthogroups": [orthogroup], "species": [spec]}
        else:
            genes[gene]["orthogroups"].append(orthogroup)
            genes[gene]["species"].append(spec)

            with open(snakemake.log[0], "a") as logfile:
                logfile.write(f'WARNING: gene {gene} already associated to orthogroup \
{genes[gene]["orthogroups"]}. \n')
            # genes[gene]["orthogroup"].append(orthogroup)
            # genes[gene]["species"].append(spec)

        if spec not in species.keys():
            species[spec] = {"orthogroups": [orthogroup], "genes": [gene]}
        else:
            species[spec]["orthogroups"].append(orthogroup)
            species[spec]["genes"].append(gene)

        if orthogroup not in orthogroups_2_species_2_genes.keys():
            orthogroups_2_species_2_genes[orthogroup] = {spec: [gene]}
        elif spec not in orthogroups_2_species_2_genes[orthogroup].keys():
            orthogroups_2_species_2_genes[orthogroup].update({spec: [gene]})
        else:
            orthogroups_2_species_2_genes[orthogroup][spec].append(gene)


# Writing out orthology table Info statistics file
n_species = len(species.keys())
n_orthogroups = len(orthogroups.keys())
n_genes = len(genes.keys())

output_file_info.write('n_species\tn_orthogroups\tn_genes\n')
output_file_info.write(f'{n_species}\t{n_orthogroups}\t{n_genes}\n')
output_file_info.write('\nRecognised species:\n')
output_file_info.write('\n'.join(species.keys()))

# Process output files
# Saving dictionaries as pickle files so that they can be easily used by other python scripts later on.
pickle.dump(orthogroups, output_file_orthogroups, protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(genes, output_file_genes, protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(species, output_file_species, protocol=pickle.HIGHEST_PROTOCOL)
pickle.dump(orthogroups_2_species_2_genes, output_file_orthogroups_2_species_2_genes, protocol=pickle.HIGHEST_PROTOCOL)

# Close files
output_file_orthogroups.close()
output_file_genes.close()
output_file_info.close()
output_file_orthogroups_2_species_2_genes.close()
