# Compute orthogroup duplicability (proportion of species with gene duplicates)

import pickle

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
output_file_orthogroups = open(snakemake.output[0], 'w')

# Process output files
output_file_orthogroups.write('orthogroup' + '\t' 'DUP' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    unique_species = sorted(set(orthogroups[orthogroup]["species"]))
    n_unique_species = len(unique_species)
    n_species_with_duplicates = 0
    for spec in unique_species:
        if orthogroups[orthogroup]["species"].count(spec) > 1:
            n_species_with_duplicates += 1
    DUP = n_species_with_duplicates / n_unique_species
    output_file_orthogroups.write(orthogroup + '\t' + str(DUP) + '\n')

# Close files
output_file_orthogroups.close()
