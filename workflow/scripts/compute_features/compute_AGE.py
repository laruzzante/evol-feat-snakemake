# Compute orthogroup age (branchlength of orthogroup's Most Recent Common Ancestor branchlength)

import pickle
from collections import defaultdict
from MRCA_functions import get_MRCA_branchlength_from_species_list

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
MRCA_branchlengths = open(snakemake.input.MRCA_branchlengths)
output_file_orthogroups = open(snakemake.output[0], 'w')


lines = MRCA_branchlengths.readlines()

MRCA_branchlengths_dict = defaultdict(defaultdict)

for line in lines[1:]:
    splitline = line.strip().split()
    spec1 = splitline[0]
    spec2 = splitline[1]
    AGE = float(splitline[2])
    MRCA_branchlengths_dict[spec1][spec2] = AGE


# Process output files
output_file_orthogroups.write('orthogroup' + '\t' + 'AGE' + '\n')
# output_file_genes.write('gene' + '\t' 'AGE' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    species_list = set(orthogroups[orthogroup]["species"])
    AGE = get_MRCA_branchlength_from_species_list(species_list, MRCA_branchlengths_dict)
    output_file_orthogroups.write(orthogroup + '\t' + str(AGE) + '\n')
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
MRCA_branchlengths.close()
output_file_orthogroups.close()
