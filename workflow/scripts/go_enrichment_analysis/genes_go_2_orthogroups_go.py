'''
    Script which converts gene-specific go terms to orthogroup-specific go terms.
    Orthogroups go terms are assigned by merging together all of the orthgroup's gene
    go terms. Alternative ways can be implemented, such as considering only go terms
    which appear e.g. at least two times across the orthogroup's genes list.
'''

import pickle

genes = pickle.load(open(snakemake.input.genes, 'rb'))
genes_go_universe = snakemake.input.genes_go_universe
output_file = snakemake.output.orthogroups_go_universe
log_file = snakemake.log[0]

orthogroups_go_dict = {}
missing_genes_dict = {}
missing_genes_counts = 0
matching_species = {}
matching_genes = {}

with open(genes_go_universe) as f, open(log_file, 'w') as logf:
    for line in f:
        gene = line.strip().split('\t')[0]
        go_term = line.strip().split('\t')[1]
        orthogroups = []
        if gene in genes.keys():
            # checking if all genes are covered in the go_universe file
            matching_genes[gene] = ''
            # checking if all species are covered in the go_universe file
            species = set(genes[gene]["species"])
            for spec in species:
                matching_species[spec] = ''
            # extracting the orthogroup(s) to which the gene is mapped.
            orthogroups = set(genes[gene]["orthogroups"])
        else:
            # CrowGO gene universe only lists one go term at a time, hence multiple lines with same gene id are present. We only count unique gene identifiers.
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

n_matching_genes = len(matching_genes.keys())
print(f"Parsed {n_matching_genes} unique matching genes from the go_universe input file.\n")
n_orthology_genes = len(genes.keys())
ratio = 100 * n_matching_genes / n_orthology_genes
n_matching_species = len(matching_species.keys())
print(f"Parsed GO terms for {n_matching_genes} genes of {n_matching_species} species out of {n_orthology_genes} total input orthologous genes, {ratio} % orthologous gene ids recovery ratio.\n")

if missing_genes_counts > 0:
    print(f"WARNING: {missing_genes_counts} gene ids from go_universe input file not present in the input orthology table. Full list written in rule's logfile.")

with open(output_file, 'w') as f:
    for orthogroup in sorted(orthogroups_go_dict.keys()):
        # Only keeping GO terms which appear at least 1 time across the orthogroup's genes list
        filtered_go_terms_list = [x for x in orthogroups_go_dict[orthogroup] if orthogroups_go_dict[orthogroup].count(x) > 1]
        line = orthogroup + '\t' + ', '.join(filtered_go_terms_list) + '\n'
        f.write(line)
