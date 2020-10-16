# Compute orthogroup universality (relative species-span)

import pickle

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input[0], 'rb'))
output_file_orthogroups = open(snakemake.output["UNI_orthogroups"], 'w')
output_file_genes = open(snakemake.output["UNI_genes"], 'w')
n_species_from_config = float(snakemake.config["n_species"])

# Process output files
for orthogroup in sorted(orthogroups.keys()):
    UNI = len(set(orthogroups[orthogroup]["species"])) / n_species_from_config
    output_file_orthogroups.write(orthogroup + '\t' + str(UNI) + '\n')
    for gene in orthogroups[orthogroup]["genes"]:
        output_file_genes.write(gene + '\t' + str(UNI) + '\n')

output_file_orthogroups.close()
output_file_genes.close()
