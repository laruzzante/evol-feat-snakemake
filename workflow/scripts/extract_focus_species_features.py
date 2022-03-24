import pickle

merged_orthogroup_features = snakemake.input.merged_orthogroup_features[0]
orthogroup_features_by_gene = snakemake.input.orthogroup_features_by_gene[0]
orthogroups_2_species_2_genes = pickle.load(open(snakemake.input.orthogroups_2_species_2_genes, 'rb'))
outfile_orthogroups = snakemake.output.spec_orthogroups[0]
outfile_genes = snakemake.output.spec_orthogroups[1]
spec = snakemake.wildcards.spec

with open(merged_orthogroup_features) as infile, open(outfile_orthogroups, 'w') as outfile:
    print(outfile)
    header = infile.readline()
    outfile.write(header)
    for line in infile:
        orthogroup = line.split('\t')[0]
        if spec in orthogroups_2_species_2_genes[orthogroup].keys():
            outfile.write(line)

with open(orthogroup_features_by_gene) as infile, open(outfile_genes, 'w') as outfile:
    header = infile.readline()
    outfile.write(header)
    for line in infile:
        gene = line.split('\t')[0]
        sp = line.split('\t')[1]
        orthogroup = line.split('\t')[2]
        if sp == spec:
            if gene in orthogroups_2_species_2_genes[orthogroup][spec]:
                outfile.write(line)
