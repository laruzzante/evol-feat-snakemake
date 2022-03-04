# Compute orthogroup maximum synteny

import pickle

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
synteny_counts = snakemake.input.synteny_counts
output_file_orthogroups = open(snakemake.output[0], 'w')


max_synteny = {}

with open(synteny_counts) as f:
    next(f)
    for line in f:
        splitline = line.strip().split('\t')
        orthogroup = splitline[0]
        max_syn = float(splitline[2])
        gff_n_species = float(splitline[1])
        # We take 1 species of gff_n_species because we want to reach the value
        # of 1 for when there is full-synteny across species. Because synteny
        # is computed as all_species against all_species, but without the self-species
        # computation (i.e. we don't consider 1 species being syntenic with itself),
        # then we have to remove 1 species from the represented species counts.
        if gff_n_species == 1.0:
            MSY = 'NA'
        else:
            MSY = max_syn / (gff_n_species - 1)
        max_synteny[orthogroup] = MSY


# Process output files

output_file_orthogroups.write('orthogroup' + '\t' + 'MSY' + '\n')

for orthogroup in sorted(orthogroups.keys()):
    if orthogroup in max_synteny.keys():
        MSY = max_synteny[orthogroup]
        output_file_orthogroups.write(orthogroup + '\t' + str(MSY) + '\n')
    else:
        print(f"WARNING: {orthogroup} has no maximum synteny (MSY) value.")


# Close files
output_file_orthogroups.close()
