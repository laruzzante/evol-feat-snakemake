# Compute orthogroup relative synteny

import pickle
from collections import defaultdict
from MRCA_functions import get_MRCA_ntips_from_species_list

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
synteny_counts = snakemake.input.synteny_counts
MRCA_ntips = snakemake.input.MRCA_ntips
output_file_orthogroups = open(snakemake.output[0], 'w')


MRCA_ntips_dict = defaultdict(defaultdict)

with open(MRCA_ntips) as f:
    next(f)
    for line in f:
        splitline = line.strip().split()
        spec1 = splitline[0]
        spec2 = splitline[1]
        ntips = float(splitline[2])
        MRCA_ntips_dict[spec1][spec2] = ntips

average_synteny = {}

with open(synteny_counts) as f:
    next(f)
    for line in f:
        splitline = line.strip().split('\t')
        orthogroup = splitline[0]
        avg_syn = float(splitline[3])
        average_synteny[orthogroup] = avg_syn


# Process output files

output_file_orthogroups.write('orthogroup' + '\t' + 'RSY' + '\n')

for orthogroup in sorted(orthogroups.keys()):
    if orthogroup in average_synteny.keys():
        species_list = set(orthogroups[orthogroup]["species"])
        # We take 1 species of gff_n_species because we want to reach the value
        # of 1 for when there is full-synteny across species. Because synteny
        # is computed as all_species against all_species, but without the self-species
        # computation (i.e. we don't consider 1 species being syntenic with itself),
        # then we have to remove 1 species from the represented species counts.
        n_species = get_MRCA_ntips_from_species_list(species_list, MRCA_ntips_dict) - 1
        RSY = average_synteny[orthogroup] / n_species
        output_file_orthogroups.write(orthogroup + '\t' + str(RSY) + '\n')
    else:
        print(f"WARNING: {orthogroup} has no relative synteny (RSY) value.")


# Close files
output_file_orthogroups.close()
