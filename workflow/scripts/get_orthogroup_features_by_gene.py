import pickle
import statistics

merged_orthogroup_features = snakemake.input.merged_orthogroup_features
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
genes = pickle.load(open(snakemake.input.orthogroups, 'rb'))
outfile = snakemake.output

with open(merged_orthogroup_features) as f, open(outfile, 'w') as f2:
    header = f.readline()
    new_header = 'gene\t\species\t' + header
    f2.write(new_header)
    next(f)
    for line in f:
        orthogroup = line.split('\t')[0]
        genes = orthogroups['genes']
        for gene in genes:
            species = genes['species']
            for spec in species:
                new_line = gene + '\t' + species + '\t' + line
                f2.write(new_line)
