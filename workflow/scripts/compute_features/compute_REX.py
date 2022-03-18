# Compute orthogroup relative expansions (number of CAFE expansion events
# divided by total number of internal nodes + leaves arising from orthogroup's
# Most Recent Common Ancestor). Or simply the number of "events" after the MRCA

import pickle
from collections import defaultdict
from MRCA_functions import get_MRCA_ntips_from_species_list

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
MRCA_ntips = open(snakemake.input.MRCA_ntips)
expansions = snakemake.input.expansions[0]
output_file_orthogroups = open(snakemake.output[0], 'w')

expansions_dict = {}
with open(expansions) as f:
    next(f)
    for line in f:
        orthogroup = line.strip().split('\t')[0]
        exp = line.strip().split('\t')[1]
        expansions_dict[orthogroup] = float(exp)

lines = MRCA_ntips.readlines()

MRCA_ntips_dict = defaultdict(defaultdict)

for line in lines[1:]:
    splitline = line.strip().split()
    spec1 = splitline[0]
    spec2 = splitline[1]
    ntips = float(splitline[2])
    MRCA_ntips_dict[spec1][spec2] = ntips


# Process output files
output_file_orthogroups.write('orthogroup' + '\t' + 'REX' + '\n')

for orthogroup in sorted(orthogroups.keys()):
    if(orthogroup in expansions_dict.keys()):
        species_list = set(orthogroups[orthogroup]["species"])
        # In a bifurcating tree, the number of Events = 2 * (n_tips - 1)
        n_nodes = 2 * ( get_MRCA_ntips_from_species_list(species_list, MRCA_ntips_dict) - 1 )
        REX = expansions_dict[orthogroup] / n_nodes
        output_file_orthogroups.write(orthogroup + '\t' + str(REX) + '\n')
    else:
        print(f'WARNING: {orthogroup} not present in expansions counts table.')


# Close files
MRCA_ntips.close()
output_file_orthogroups.close()
