import pickle

genes = pickle.load(open(snakemake.input.genes, 'rb'))
genes_go_universe = snakemake.input.genes_go_universe
output_file = open(snakemake.output.orthogroups_go_universe, 'w')

orthogroups_go_dict = {}

with open(genes_go_universe) as f:
    for line in f:
        gene = line.strip().split('\t')[0]
        go_term = line.strip().split('\t')[1]
        if gene in genes.keys():
            orthogroups = genes[gene]["orthogroups"]
        else:
            print(f"WARNING: geneid '{gene}' from go_universe not present in orthology table. Skipping.")
            next(f)
        for orthogroup in orthogroups:
            if orthogroup in orthogroups_go_dict.keys():
                orthogroups_go_dict[orthogroup].append(go_term)
            else:
                orthogroups_go_dict[orthogroup] = [go_term]

with open(output_file, 'w') as f:
    for orthogroup in orthogroups_go_dict.keys():
        line = orthogroup + '\t' + ', '.join(orthogroups_go_dict[orthogroup]) + '\n'
        f.write(line)
