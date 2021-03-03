import pandas as pd

orthogroup_sets = open(snakemake.input.ogsets)
features = open(snakemake.input.features)
outfile = open(snakemake.output[0], 'w')

features_lines = features.readlines()
features_names = features_lines[0].strip().split('\t')[1:]

orthogroup_features = {}
for line in features_lines[1:]:
    split_line = line.strip().split('\t')
    orthogroup = split_line[0]
    orthogroup_features[orthogroup] = split_line[1:]

outfile.write('orthogroup' + '\t' + 'set' + '\t' + '\t'.join(features_names) + '\n')
for line in orthogroup_sets:
    split_line = line.strip().split('\t')
    setname = split_line[0]
    orthogroups = split_line[1:]
    for orthogroup in orthogroups:
        if orthogroup in orthogroup_features.keys():
            values = "\t".join(orthogroup_features[orthogroup])
            outfile.write(orthogroup + '\t' + setname + '\t' + values + '\n')
        else:
            print(f'WARNING: {orthogroup} has no feature values.')

# Close files
orthogroup_sets.close()
features.close()
outfile.close()
