import pickle

genes = pickle.load(open(snakemake.input.genes, 'rb'))
genes_go_universe = snakemake.input.genes_go_universe
output_file = snakemake.output.orthogroups_go_universe
log_file = snakemake.log[0]

orthogroups_go_dict = {}
missing_genes_dict = {}
missing_genes_counts = 0

with open(genes_go_universe) as f, open(log_file, 'w') as logf:
    for line in f:
        gene = line.strip().split('\t')[0]
        go_term = line.strip().split('\t')[1]
        if gene in genes.keys():
            orthogroups = genes[gene]["orthogroups"]
        else:
            # CrowGO gene universe lists only one goterm at a time, hence multiple lines with same gene id are present. We only count unique gene identifiers.
            if gene not in missing_genes_dict.keys():
                missing_genes_dict[gene] = ''
                missing_genes_counts += 1
                logf.write(f"WARNING: geneid '{gene}' from the go_universe input file are not present in orthology table. Skipping.\n")
            next(f)
        for orthogroup in orthogroups:
            if orthogroup in orthogroups_go_dict.keys():
                orthogroups_go_dict[orthogroup].append(go_term)
            else:
                orthogroups_go_dict[orthogroup] = [go_term]

if missing_genes_counts > 0:
    print(f"WARNING: {missing_genes_counts} gene ids from go_universe file not present in orthology table. Full list written in rule's logfile.")

with open(output_file, 'w') as f:
    for orthogroup in sorted(orthogroups_go_dict.keys()):
        line = orthogroup + '\t' + ', '.join(orthogroups_go_dict[orthogroup]) + '\n'
        f.write(line)
