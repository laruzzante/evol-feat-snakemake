import pickle

orthogroupsets_file = open(snakemake.input[0])
orthogroups_ortho_dict = pickle.load(open(snakemake.input.orthogroups, 'rb'))
output_file = open(snakemake.output[0], 'w')

ogsets = {}

for line in orthogroupsets_file:
    split_line = line.strip().split('\t')
    orthogroup_id = split_line[0]
    set_name = split_line[1]
    if orthogroup_id in orthogroups_ortho_dict.keys():
        if set_name in ogsets.keys():
            ogsets[set_name].append(orthogroup_id)
        else:
            ogsets[set_name] = [orthogroup_id]
    else:
        print(f'WARNING: orthogroup {orthogroup_id} missing from orthology table. Discarding.')

for set in ogsets.keys():
    output_file.write(set + '\t' + '\t'.join(ogsets[set]) + '\n')

# Close files
orthogroupsets_file.close()
output_file.close()
