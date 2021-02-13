# Compute orthogroup copy-number variation (per species standard deviation of gene copies, divided by average copy-number)

import pickle
import statistics as stat

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input[0], 'rb'))
ACN_orthogroups = open(snakemake.input["ACN_orthogroups"])
output_file_orthogroups = open(snakemake.output["CNV_orthogroups"], 'w')
output_file_genes = open(snakemake.output["CNV_genes"], 'w')

# Create ACN orthogroups dictionary from output of compute_average_copy_number rule:
ACN_dictionary = {}
next(ACN_orthogroups) # skipping first line as it is the file header
for line in ACN_orthogroups:
    orthogroup = line.split('\t')[0]
    ACN = line.split('\t')[1]
    ACN_dictionary[orthogroup] = float(ACN)

# Process output files
output_file_orthogroups.write('orthogroup' + '\t' 'CNV' + '\n')
output_file_genes.write('gene' + '\t' 'CNV' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    unique_species = sorted(set(orthogroups[orthogroup]["species"]))
    species_copy_counts = []
    for spec in unique_species:
        species_copy_count = orthogroups[orthogroup]["species"].count(spec)
        species_copy_counts.append(species_copy_count)
    species_copy_counts_stdev = stat.stdev(species_copy_counts)
    CNV = species_copy_counts_stdev / ACN_dictionary[orthogroup]
    output_file_orthogroups.write(orthogroup + '\t' + str(CNV) + '\n')
    for gene in sorted(orthogroups[orthogroup]["genes"]):
        output_file_genes.write(gene + '\t' + str(CNV) + '\n')

# Close files
output_file_orthogroups.close()
output_file_genes.close()
