# Compute orthogroup duplicability (proportion of species with gene duplicates)

import pickle

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input[0], 'rb'))
output_file_orthogroups = open(snakemake.output["DUP_orthogroups"], 'w')
output_file_genes = open(snakemake.output["DUP_genes"], 'w')

# Process output files
output_file_orthogroups.write('orthogroup' + '\t' 'DUP' + '\n')
output_file_genes.write('gene' + '\t' 'DUP' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    unique_species = sorted(set(orthogroups[orthogroup]["species"]))
    n_unique_species = len(unique_species)
    n_species_with_duplicates = 0
    for spec in unique_species:
        if orthogroups[orthogroup]["species"].count(spec) > 1:
            n_species_with_duplicates += 1
    DUP = n_species_with_duplicates / n_unique_species
    output_file_orthogroups.write(orthogroup + '\t' + str(DUP) + '\n')
    for gene in sorted(orthogroups[orthogroup]["genes"]):
        output_file_genes.write(gene + '\t' + str(DUP) + '\n')

# Close files
output_file_orthogroups.close()
output_file_genes.close()
