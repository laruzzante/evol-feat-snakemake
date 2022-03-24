import pickle
import statistics

merged_orthogroup_features = snakemake.input.merged_orthogroup_features[0]
orthogroup_features_by_gene = snakemake.input.orthogroup_features_by_gene[0]
species = pickle.load(open(snakemake.input.species, 'rb'))
outfile_orthogroups = snakemake.output.spec_orthogroups[0]
outfile_genes = snakemake.output.spec_orthogroups[0]
spec = snakemake.params.spec

with open(merged_orthogroup_features) as f, open(outfile_orthogroups, 'w') as f2:
    header = f.readline()
    f2.write(header)
    for line in f:
        orthogroup = line.split('\t')[0]
        if orthogroup in species[spec]['orthogroups']:
            f2.write(line)

with open(orthogroup_features_by_gene) as f, open(outfile_genes, 'w') as f2:
    header = f.readline()
    f2.write(header)
    for line in f:
        gene = line.split('\t')[0]
        if orthogroup in species[spec]['genes']:
            f2.write(line)
