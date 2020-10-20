# Compute orthogroup copy-number variation (per species standard deviation of gene copies, divided by average copy-number)

import pickle
import statistics as stat

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input[0], 'rb'))
output_file_orthogroups = open(snakemake.output["CNV_orthogroups"], 'w')
output_file_genes = open(snakemake.output["CNV_genes"], 'w')

# Process output files
output_file_orthogroups.write('orthogroup' + '\t' 'CNV' + '\n')
output_file_genes.write('gene' + '\t' 'CNV' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    unique_species = sorted(set(orthogroups[orthogroup]["species"]))
    n_unique_species = len(unique_species)
    n_genes = len(set(orthogroups[orthogroup]["genes"]))
    ACN = n_genes / n_unique_species # average copy-number
    species_copy_counts = []
    for spec in unique_species:
        species_copy_count = orthogroups[orthogroup]["species"].count(spec)
        species_copy_counts.append(species_copy_count)
    species_copy_counts_stdev = stat.stdev(species_copy_counts)
    CNV = species_copy_counts_stdev / ACN
    output_file_orthogroups.write(orthogroup + '\t' + str(CNV) + '\n')
    for gene in sorted(orthogroups[orthogroup]["genes"]):
        output_file_genes.write(gene + '\t' + str(CNV) + '\n')

# Close files
output_file_orthogroups.close()
output_file_genes.close()
