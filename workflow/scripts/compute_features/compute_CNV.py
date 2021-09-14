# Compute orthogroup copy-number variation (per species standard deviation of gene copies, divided by average copy-number)

import pickle
import statistics as stat

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
# ACN_orthogroups = open(snakemake.input.ACN_orthogroups[0])
output_file_orthogroups = open(snakemake.output[0], 'w')

# Create ACN orthogroups dictionary from output of compute_average_copy_number rule:
# ACN_dictionary = {}
# next(ACN_orthogroups) # skipping first line as it is the file header
# for line in ACN_orthogroups:
#     orthogroup = line.split('\t')[0]
#     ACN = line.split('\t')[1]
#     ACN_dictionary[orthogroup] = float(ACN)
#
# # Process output files
# output_file_orthogroups.write('orthogroup' + '\t' 'CNV' + '\n')
# for orthogroup in sorted(orthogroups.keys()):
#     unique_species = sorted(set(orthogroups[orthogroup]["species"]))
#     species_copy_counts = []
#     for spec in unique_species:
#         species_copy_count = orthogroups[orthogroup]["species"].count(spec)
#         species_copy_counts.append(species_copy_count)
#     species_copy_counts_stdev = stat.stdev(species_copy_counts)
#     CNV = species_copy_counts_stdev / ACN_dictionary[orthogroup]
#     output_file_orthogroups.write(orthogroup + '\t' + str(CNV) + '\n')


#######################


ACN_dictionary = {}
output_file_orthogroups.write('orthogroup' + '\t' 'CNV' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    # Computing average copy-number first
    unique_species = sorted(set(orthogroups[orthogroup]["species"]))
    n_unique_species = len(unique_species)
    n_genes = len(set(orthogroups[orthogroup]["genes"]))
    ACN = n_genes / n_unique_species
    ACN_dictionary[orthogroup] = ACN

    # Computing copy-number variation
    species_copy_counts = []
    for spec in unique_species:
        species_copy_count = orthogroups[orthogroup]["species"].count(spec)
        species_copy_counts.append(species_copy_count)
    if len(species_copy_counts) == 1:
        print(f'WARNING: orthogroup {orthogroup} is present in only 1 species, {spec}')
        species_copy_counts_stdev = 0
    else:
        species_copy_counts_stdev = stat.stdev(species_copy_counts)
    CNV = species_copy_counts_stdev / ACN_dictionary[orthogroup]
    output_file_orthogroups.write(orthogroup + '\t' + str(CNV) + '\n')

# Close files
output_file_orthogroups.close()
