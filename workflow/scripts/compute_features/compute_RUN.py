# Compute orthogroup age (branchlength of orthogroup's Most Recent Common Ancestor branchlength)

import pickle
from collections import defaultdict

# Retrieve information from Snakemake
orthogroups = pickle.load(open(snakemake.input[0], 'rb'))
mrca_ntips = open(snakemake.input.mrca_ntips)
output_file_orthogroups = open(snakemake.output[0], 'w')


def get_MRCA_ntips_from_specieslist(species_list, mrca_ntips_dict):

    ntips_max = 0

    for spec1 in sorted(species_list):
        for spec2 in sorted(species_list):
            if(spec1 != spec2):
                if mrca_ntips_dict[spec1][spec2]:
                    ntips = mrca_ntips_dict[spec1][spec2]
                elif mrca_ntips_dict[spec2][spec1]:
                    ntips = mrca_ntips_dict[spec2][spec1]
                else:
                    print('ERROR: species combination not present in MRCA ntips dictionary.')
                    sys.close()

                if ntips > ntips_max:
                    ntips_max = ntips

    if ntips_max == 0:
        print('ERROR: species list ' + ' '.join(species_list) + ' returns a MRCA ntips of 0.')

    return ntips_max


lines = mrca_ntips.readlines()

mrca_ntips_dict = defaultdict(defaultdict)

for line in lines[1:]:
    splitline = line.strip().split()
    spec1 = splitline[0]
    spec2 = splitline[1]
    ntips = float(splitline[2])
    mrca_ntips_dict[spec1][spec2] = ntips


# Process output files
output_file_orthogroups.write('orthogroup' + '\t' 'RUN' + '\n')
# output_file_genes.write('gene' + '\t' 'AGE' + '\n')
for orthogroup in sorted(orthogroups.keys()):
    species_list = set(orthogroups[orthogroup]["species"])
    RUN = len(species_list) / get_MRCA_ntips_from_specieslist(species_list, mrca_ntips_dict)
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
mrca_ntips.close()
output_file_orthogroups.close()
