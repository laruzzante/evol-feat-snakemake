# Compute orthogroup age (branchlength of orthogroup's Most Recent Common Ancestor branchlength)

import pickle
from collections import defaultdict

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input[0], 'rb'))
mrca_branchlengths = open(snakemake.input.mrca_Ntip)
output_file_orthogroups = open(snakemake.output[0], 'w')


def getRUN_from_specieslist(species_list, AGE_LCAbranchlengths_dict):

    AGE_max = 0

    for spec1 in sorted(species_list):
        for spec2 in sorted(species_list):
            if(spec1 != spec2):
                if AGE_LCAbranchlengths_dict[spec1][spec2]:
                    AGE = AGE_LCAbranchlengths_dict[spec1][spec2]
                elif AGE_LCAbranchlengths_dict[spec2][spec1]:
                    AGE = AGE_LCAbranchlengths_dict[spec2][spec1]
                else:
                    print('ERROR: species combination not present in LCA branchlengths dictionary.')
                    sys.close()

                if AGE > AGE_max:
                    AGE_max = AGE

    if AGE_max == 0:
        print('ERROR: species list ' + ' '.join(species_list) + ' returns a LCA branchlengths of 0.')

    return AGE_max


lines = mrca_branchlengths.readlines()

AGE_LCAbranchlengths_dict = defaultdict(defaultdict)

for line in lines[1:]:
    splitline = line.strip().split()
    spec1 = splitline[0]
    spec2 = splitline[1]
    AGE = float(splitline[2])
    AGE_LCAbranchlengths_dict[spec1][spec2] = AGE


# Process output files
output_file_orthogroups.write('orthogroup' + '\t' 'AGE' + '\n')
# output_file_genes.write('gene' + '\t' 'AGE' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    species_list = set(orthogroups[orthogroup]["species"])
    AGE = getAGE_from_specieslist(species_list, AGE_LCAbranchlengths_dict)
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
mrca_branchlengths.close()
output_file_orthogroups.close()
