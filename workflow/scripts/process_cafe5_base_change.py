# Create table of orthogroup copy number variations across tree from the CAFE5 Base_change results table
# Each line of the CAFE5 Base_change results table contains the ortogroup id in the first column and
# subsequently its copy-number variation count at each tree node.


# Retrieve information from Snakemake
input_file = snakemake.input.cafe_results
output_file = open(snakemake.output.copy_number_variation_table, 'w')

# Initiliazing dictionary where we will store CAFE events, i.e.:
# Counts of Expansions, Stabilities or Contractions
copy_number_variations = {}

with open(input_file) as f:
    next(f) # skipping header row of cafe5 Base_change table

    for line in f:

        n_expansions = 0
        n_stabilities = 0
        n_contractions = 0

        splitted_line = line.strip().split('\t')
        orthogroup = splitted_line[0]

        for element in splitted_line[1:]:
            n = int(element)
            if n > 0:
                n_expansions += 1
            elif n == 0:
                n_stabilities += 1
            elif n < 0:
                n_contractions += 1
            else:
                print(f"WARNING: {orthogroup} has no valid copy-number integer value in the cafe results table.")

        copy_number_variations[orthogroup] = {'EXP': n_expansions,
                                              'STA': n_stabilities,
                                              'CON': n_contractions}


# Process output files
output_file.write('orthogroup\tEXP\tSTA\tCON\n')

for orthogroup in sorted(copy_number_variations.keys()):
    output_file.write(f"{orthogroup}\t{copy_number_variations[orthogroup]['EXP']}\t{copy_number_variations[orthogroup]['STA']}\t{copy_number_variations[orthogroup]['CON']}\n")

# Close files
output_file.close()
