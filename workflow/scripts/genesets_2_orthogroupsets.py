import pickle

genesets_file = open(snakemake.input[0])
genes_ortho_dict = pickle.load(open(snakemake.input.genes, 'rb'))
output_file = open(snakemake.output[0], 'w')

genesets = {}
ogsets = {}

for line in genesets_file:
    split_line = line.strip().split('\t')
    set_name = split_line[0]
    genes = split_line[1:]
    genesets[set_name] = genes

for set in genesets.keys():
    for gene in genesets[set]:
        if gene in genes_ortho_dict.keys():
            og = genes_ortho_dict[gene]['orthogroup']
            if set not in ogsets.keys():
                ogsets[set] = [og]
            else:
                ogsets[set].append(og)
        else:
            print(f'WARNING: gene {gene} missing from orthology table. Discarding.')

for set in ogsets.keys():
    output_file.write(set + '\t' + '\t'.join(ogsets[set]) + '\n')

# Close files
genesets_file.close()
output_file.close()
