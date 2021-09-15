# Compute orthogroup age (branchlength of orthogroup's Most Recent Common Ancestor branchlength)

import pickle
from collections import defaultdict
from MRCA_functions import get_MRCA_ntips_from_species_list

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
MRCA_ntips = open(snakemake.input.MRCA_ntips)
output_file_orthogroups = open(snakemake.output[0], 'w')


lines = MRCA_ntips.readlines()

MRCA_ntips_dict = defaultdict(defaultdict)

for line in lines[1:]:
    splitline = line.strip().split()
    spec1 = splitline[0]
    spec2 = splitline[1]
    ntips = float(splitline[2])
    MRCA_ntips_dict[spec1][spec2] = ntips


# Process output files
output_file_orthogroups.write('orthogroup' + '\t' + 'RUN' + '\n')
# output_file_genes.write('gene' + '\t' 'AGE' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    species_list = set(orthogroups[orthogroup]["species"])
    RUN = len(species_list) / get_MRCA_ntips_from_species_list(species_list, MRCA_ntips_dict)
    output_file_orthogroups.write(orthogroup + '\t' + str(RUN) + '\n')
    # for gene in sorted(orthogroups[orthogroup]["genes"]):
    #     output_file_genes.write(gene + '\t' + str(AGE) + '\n')

        # Counter idea, but computation is fast anyways, so optimization not so useful:
        # i += 1
        # if (i / len(ogDict) > len(ogDict) * 0.25) and (i / len(ogDict) <= len(ogDict) * 0.25 + 1):
        #     print('\n\t... 25%')
        # elif (i / len(ogDict) > len(ogDict) * 0.5) and (i / len(ogDict) <= len(ogDict) * 0.5 + 1):
        #     print('\t... 50%')
        # elif (i / len(ogDict) > len(ogDict) * 0.75) and (i / len(ogDict) <= len(ogDict) * 0.75 + 1):
        #     print('\t... 75%')
        # elif i == len(ogDict):
        #     print('\t... 100%')

# Close files
MRCA_ntips.close()
output_file_orthogroups.close()
