# Retrieve information from Snakemake
input_file = open(snakemake.input[0])
output_file = open(snakemake.output[0], 'w')
n_lines = len(input_file)

# Process file
output_file.write(n_lines)
