import pickle
import statistics

merged_orthogroup_features = snakemake.input.merged_orthogroup_features[0]
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
genes = pickle.load(open(snakemake.input.genes, 'rb'))
outfile = snakemake.output[0]

with open(merged_orthogroup_features) as f, open(outfile, 'w') as f2:
    header = f.readline()
    new_header = 'gene\tspecies\t' + header
    f2.write(new_header)
    for line in f:
        orthogroup = line.split('\t')[0]
        genes_list = orthogroups[orthogroup]['genes']
        for gene in genes_list:
            species = genes[gene]['species']
            for spec in species:
                new_line = gene + '\t' + spec + '\t' + line
                f2.write(new_line)
