import pickle

genes = pickle.load(open(snakemake.input.genes, 'rb'))
concatenated_gff_file = snakemake.input.gff[0]
output_ordered_gff_genes = open(snakemake.output.ordered_gff_genes, 'w')

num_lines = sum(1 for line in open(concatenated_gff_file))
i = 0

gene_coordinates = {}
gene_centers = {}
# ordered_gff_genes = {}

with open(concatenated_gff_file) as infile:
    print('Processing concatenated GFF file: ' + concatenated_gff_file)
    for line in infile:

        i += 1
        if i%50000 == 0:
            progress = round(100*i/num_lines , 2)
            print('\t\t' + str(progress) + ' %')

        splitline = line.strip().split('\t')
        geneid = splitline[9]
        spec = splitline[10]
        gene_start = int(splitline[3])
        gene_end = int(splitline[4])

        if spec not in gene_coordinates.keys():
            gene_coordinates[spec] = {}
        if geneid not in gene_coordinates[spec].keys():
            # Not using defaultdict because they are not easily handled by Pickle
            # Hence, the update method.
            gene_coordinates[spec].update({geneid: {'start': gene_start, 'end': gene_end}})
        else:
            if gene_start < gene_coordinates[spec][geneid]['start']:
                gene_coordinates[spec][geneid]['start'] = gene_start
            if gene_end > gene_coordinates[spec][geneid]['end']:
                gene_coordinates[spec][geneid]['end'] = gene_end

for spec in gene_coordinates.keys():
    for geneid in gene_coordinates[spec].keys():
        gene_start = gene_coordinates[spec][geneid]['start']
        gene_end = gene_coordinates[spec][geneid]['end']
        gene_center = 0.5 * (gene_start + gene_end)
        if spec not in gene_centers.keys():
            gene_centers[spec] = {}
        if gene_center not in gene_centers[spec].keys():
            # List because some gene centers may not be unique to only one gene in the GFF file.
            gene_centers[spec].update({gene_center: [geneid]})
        else:
            gene_centers[spec][gene_center].append(geneid)

output_ordered_gff_genes.write('spec\torthogroup\tgene\tgene_center\n')
for spec in sorted(gene_centers.keys()):
    for gene_center in sorted(gene_centers[spec].keys()):
        for geneid in sorted(gene_centers[spec][gene_center]):
            ogs = list(set(genes[geneid]["orthogroups"]))
            # Ideally there should be only 1 OG per gene, however because OrthoDB
            # recently changed this to take care of cut genes mapping to several OGs,
            # we have to loop through all OGs. Ideally this should not be the case.
            # This applies to 2530 genes in the 170 arthropoda species odbv10.1 dataset.
            for og in ogs:
                output_ordered_gff_genes.write(spec + '\t' + og + '\t' + geneid + '\t' + str(gene_center) + '\n')
                # if spec not in ordered_gff_genes.keys():
                #     ordered_gff_genes[spec] = [geneid]
                # else:
                #     ordered_gff_genes[spec].append(geneid)

# Close files
output_ordered_gff_genes.close()
