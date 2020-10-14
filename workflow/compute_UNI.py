# Computing orthogroup universality

# Retrieve information from Snakemake
input_file = open(snakemake.input[0])
output_file = open(snakemake.output[0], 'w')



# Process file
output_file.write(input_file.readline())
