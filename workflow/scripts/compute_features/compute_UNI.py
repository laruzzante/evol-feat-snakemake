# Compute orthogroup universality (relative species-span)

import pickle

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
info =  open(snakemake.input.info)
output_file_orthogroups = open(snakemake.output[0], 'w')
# output_file_genes = open(snakemake.output["UNI_genes"], 'w')

info_lines = info.readlines()
n_species = int(info_lines[1].split('\t')[0])

# n_species_from_config = float(snakemake.config["n_species"])

# Process output files
output_file_orthogroups.write('orthogroup' + '\t' + 'UNI' + '\n')
# output_file_genes.write('gene' + '\t' 'UNI' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    n_unique_species = len(set(orthogroups[orthogroup]["species"]))
    UNI =  n_unique_species / n_species
    output_file_orthogroups.write(orthogroup + '\t' + str(UNI) + '\n')
    # for gene in sorted(orthogroups[orthogroup]["genes"]):
    #     output_file_genes.write(gene + '\t' + str(UNI) + '\n')

# Close files
info.close()
output_file_orthogroups.close()
# output_file_genes.close()
