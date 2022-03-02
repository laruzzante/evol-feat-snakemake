import pickle

genes = pickle.load(open(snakemake.input.genes, 'rb'))
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
concatenated_gff_file = snakemake.input.gff[0]
output_file = open(snakemake.output.allspecies_synteny_scores, 'w')


num_lines = sum(1 for line in open(concatenated_gff_file))
i = 0

gene_coordinates = {}
gene_centers = {}
ordered_gff_genes = {}

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

for spec in gene_centers.keys():
    for gene_center in sorted(gene_centers[spec].keys()):
        for geneid in gene_centers[spec][gene_center]:
            if spec not in ordered_gff_genes.keys():
                ordered_gff_genes[spec] = [geneid]
            else:
                ordered_gff_genes[spec].append(geneid)


synteny_dict = {}

# Looping through all of the species so to compute SYN for each species as a reference species and the final
# SYN score will be the average of all.
k = 0
for ref_spec in ordered_gff_genes.keys():

    k += 1
    n_ref_spec = len(ordered_gff_genes.keys())
    print('Computing SYN with ' + ref_spec + ' as reference species. ' + str(k) + '/' + str(n_ref_spec))

    i = 0
    for ref_geneid in ordered_gff_genes[ref_spec]:
        n_genes = len(ordered_gff_genes[ref_spec])
        if i%100 == 0:
            progress = round(100*i/n_genes, 2)
            print('... ' + str(progress) + ' %')
        if ref_geneid in genes.keys():
             # Ideally there should be only 1 OG per gene, however because OrthoDB
             # recently changed this to take care of cut genes mapping to several OGs,
             # we have to loop through all OGs. Ideally this should not be the case.
            for og in genes[ref_geneid]["orthogroups"]:
                # Taking care of 1st gene
                if i == 0:
                    ref_geneid_left = None
                else:
                    ref_geneid_left = ordered_gff_genes[ref_spec][i - 1]

                # Taking care of last gene
                if i == len(ordered_gff_genes[ref_spec]) - 1:
                    ref_geneid_right = None
                else:
                    ref_geneid_right = ordered_gff_genes[ref_spec][i + 1]

                # Finding corresponding left and right OGs
                is_og_left = False
                is_og_right = False

                if ref_geneid_left in genes.keys() and ref_geneid_right in genes.keys():
                    ogs_left = genes[ref_geneid_left]["orthogroups"]
                    ogs_right = genes[ref_geneid_right]["orthogroups"]
                    is_og_left = True
                    is_og_right = True
                elif ref_geneid_left in genes.keys() and ref_geneid_right not in genes.keys():
                    ogs_left = genes[ref_geneid_left]["orthogroups"]
                    ogs_right = None
                    is_og_left = True
                elif ref_geneid_left not in genes.keys() and ref_geneid_right in genes.keys():
                    ogs_right = genes[ref_geneid_right]["orthogroups"]
                    ogs_left = None
                    is_og_right = True
                else:
                    i += 1
                    continue

                # Checking if left and right genes also map to same left and right OGs in other species.
                SYN = 0
                syn_left = False
                syn_right = False
                og_species = sorted(orthogroups[og]["species"])
                og_species.remove(ref_spec)

                for spec in og_species:
                    geneids = orthogroups[og]["genes"]

                    for geneid in geneids:
                        # Remove genes that do not belong to the current species
                        if spec not in genes[geneid]["species"]:
                            geneids.remove(geneid)

                    if is_og_left == True:
                        for og_left in ogs_left:
                            if syn_left == True:
                                break
                            geneids_left = []
                            if spec in orthogroups[og_left]["species"] and spec in ordered_gff_genes.keys():
                                geneids_left.extend(orthogroups[og_left]["genes"])
                                for geneid in geneids:
                                    if geneid in ordered_gff_genes[spec]:
                                        geneid_index = ordered_gff_genes[spec].index(geneid)
                                        if geneid_index != 0:  # Removing the case where there is no left gene in the species (i.e. first gene)
                                            geneid_left = ordered_gff_genes[spec][geneid_index - 1]
                                            if geneid_left in geneids_left:
                                                syn_left = True
                                                break
                                        else:
                                            continue
                                    # else:
                                    #     print("WARNING: orthologous gene_id " + geneid + " not in " + spec + " GFF file!")
                                    #     continue

                    if is_og_right == True:
                        for og_right in ogs_right:
                            if syn_right == True:
                                break
                            geneids_right = []
                            if spec in orthogroups[og_right]["species"] and spec in ordered_gff_genes.keys():
                                for og_right in ogs_right:
                                    geneids_right.extend(orthogroups[og_right]["genes"])
                                    for geneid in geneids:
                                        if geneid in ordered_gff_genes[spec]:
                                            geneid_index = ordered_gff_genes[spec].index(geneid)
                                            if geneid_index < len(ordered_gff_genes[spec]) - 1:  # Removing the case where there is no right gene in the species (i.e. last gene)
                                                geneid_right = ordered_gff_genes[spec][geneid_index + 1]
                                                if geneid_right in geneids_right:
                                                    syn_right = True
                                                    break
                                            else:
                                                continue
                                        # else:
                                        #     print("WARNING: orthologous gene_id " + geneid + " not in " + spec + " GFF file!")
                                        #     continue

                    if syn_left == True or syn_right == True:
                        SYN += 1

            # Dividing SYN score by the number of species without the reference species,
            # because if OG has complete synteny with all other species, then SYN score will be 1 (in a scale from 0 to 1)
            SYN = SYN / len(og_species)
            if og not in synteny_dict.keys():
                synteny_dict[og] = {}
            if ref_spec not in synteny_dict[og].keys():
                synteny_dict[og].update({ref_spec: SYN})
            else:
                if SYN > synteny_dict[og][ref_spec]:
                    synteny_dict[og][ref_spec] = SYN

            i += 1

        else:
            i += 1
            continue

# Process output files
output_file.write('orthogroup\treference_species\tsynteny_score\n')
for orthogroup in sorted(synteny_dict.keys()):
    for ref_spec in sorted(synteny_dict[orthogroup].keys()):
        output_file.write(f"{orthogroup}\t{ref_spec}\t{synteny_dict[orthogroup][ref_spec]}\n")

# Close files
output_file.close()
