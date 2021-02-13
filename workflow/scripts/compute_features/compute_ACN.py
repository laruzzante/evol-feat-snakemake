# Compute orthogroup average copy-number (average of per species gene copy-number)

import pickle

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input[0], 'rb'))
output_file_orthogroups = open(snakemake.output["ACN_orthogroups"], 'w')
output_file_genes = open(snakemake.output["ACN_genes"], 'w')

# Process output files
output_file_orthogroups.write('orthogroup' + '\t' 'ACN' + '\n')
output_file_genes.write('gene' + '\t' 'ACN' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    unique_species = sorted(set(orthogroups[orthogroup]["species"]))
    n_unique_species = len(unique_species)
    n_genes = len(set(orthogroups[orthogroup]["genes"]))
    ACN = n_genes / n_unique_species
    output_file_orthogroups.write(orthogroup + '\t' + str(ACN) + '\n')
    for gene in sorted(orthogroups[orthogroup]["genes"]):
        output_file_genes.write(gene + '\t' + str(ACN) + '\n')

# Close files
output_file_orthogroups.close()
output_file_genes.close()
