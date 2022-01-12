# Merge all computed features and print them in tsv file

import pandas as pd
from functools import reduce

# Retrieve information from Snakemake
orthogroup_features_files = snakemake.input
output_file = open(snakemake.output[0], 'w')

# Process output files
orthogroup_features = []
for input_file in sorted(orthogroup_features_files):
    # I <open(input_file)> and don't just use <input_file> here because if given a string,
    # pandas.read_csv() will automatically close the file, which is instead needed to be open for later.
    df = pd.read_csv(open(input_file), sep='\t')
    orthogroup_features.append(df)

merged_orthogroup_features = reduce(lambda left, right: pd.merge(left, right, on='orthogroup', how='outer'), orthogroup_features).fillna('NA')
merged_orthogroup_features.to_csv(output_file, sep='\t', index=False)

# Close files
output_file.close()
