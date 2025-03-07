# Compute orthogroup average synteny

import pickle

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
synteny_counts = snakemake.input.synteny_counts
output_file_orthogroups = open(snakemake.output[0], 'w')


synteny = {}

with open(synteny_counts) as f:
    next(f)
    for line in f:
        splitline = line.strip().split('\t')
        orthogroup = splitline[0]
        avg_syn = float(splitline[3])
        gff_n_species = float(splitline[1])
        # We take 1 species of gff_n_species because we want to reach the value
        # of 1 for when there is full-synteny across species. Because synteny
        # is computed as all_species against all_species, but without the self-species
        # computation (i.e. we don't consider 1 species being syntenic with itself),
        # then we have to remove 1 species from the represented species counts.
        if gff_n_species == 1.0:
            SYN = 'NA'
        else:
            SYN = avg_syn / (gff_n_species - 1)
        synteny[orthogroup] = SYN


# Process output files

output_file_orthogroups.write('orthogroup' + '\t' + 'SYN' + '\n')

for orthogroup in sorted(orthogroups.keys()):
    if orthogroup in synteny.keys():
        SYN = synteny[orthogroup]
        output_file_orthogroups.write(orthogroup + '\t' + str(SYN) + '\n')
    else:
        print(f"WARNING: {orthogroup} has no average synteny (SYN) value.")


# Close files
output_file_orthogroups.close()
