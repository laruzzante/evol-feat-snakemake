# Compute orthogroup relative universality (species span divided by the total number of species in the phylogeny originating from orthogroup's Most Recent Common Ancestor)

import pickle
from collections import defaultdict
from MRCA_functions import get_MRCA_ntips_from_species_list

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
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


# Process output files
output_file_orthogroups.write('orthogroup' + '\t' + 'RUN' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    species_list = set(orthogroups[orthogroup]["species"])
    RUN = len(species_list) / get_MRCA_ntips_from_species_list(species_list, MRCA_ntips_dict)
    output_file_orthogroups.write(orthogroup + '\t' + str(RUN) + '\n')


# Close files
output_file_orthogroups.close()
