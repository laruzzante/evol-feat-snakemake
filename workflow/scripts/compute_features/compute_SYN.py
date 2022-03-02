# Compute orthogroup synteny

import pickle

# Retrieve information from Snakemake
genes = pickle.load(open(snakemake.input.genes, 'rb'))
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
ordered_gff_genes = pickle.load(open(snakemake.input.ordered_gff_genes, 'rb'))
output_file_orthogroups = open(snakemake.output[0], 'w')


# import statistics

synteny_dict = {} #defaultdict(lambda: defaultdict(list))

# Looping through all of the species so to compute SYN for each species as a reference species and the final
# SYN score will be the average of all.
k = 0
for ref_spec in ordered_gff_genes.keys():

    k += 1
    n_ref_spec = len(ordered_gff_genes.keys())
    print('Computing SYN with ' + ref_spec + ' as reference species. ' + str(k) + '/' + str(n_ref_spec))

    i = 0
    for ref_geneid in ordered_gff_genes[ref_spec]:
        if ref_geneid in genes.keys():
            # Finding corresponding OG
            og = genes[ref_geneid]

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
                og_left = genes[ref_geneid_left]
                og_right = genes[ref_geneid_right]
                is_og_left = True
                is_og_right = True
            elif ref_geneid_left in genes.keys() and ref_geneid_right not in genes.keys():
                og_left = genes[ref_geneid_left]
                is_og_left = True
            elif ref_geneid_left not in genes.keys() and ref_geneid_right in genes.keys():
                og_right = genes[ref_geneid_right]
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

                if is_og_left is True and spec in orthogroups[og_left]["species"]:
                    geneids_left = orthogroups[og_left]["genes"]
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

                if is_og_right is True and spec in orthogroups[og_right]["species"]:
                    geneids_right = orthogroups[og_right]["genes"]
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

                if syn_left is True or syn_right is True:
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



for og in synteny_dict.keys():
    for

for og in orthoGroupsMetrics.keys():

    if og in syntenyDict.keys():

        max_synteny_for_all_refSpecies = []

        for refSpec in syntenyDict[og].keys():
            # It can happen that for multicopy genes the og-synteny will be computed several times,
            # hence only the maximal synteny score per OG is retained.
            max_synteny_per_OG = max(syntenyDict[og][refSpec])
            max_synteny_for_all_refSpecies.append(max_synteny_per_OG)
        # These maximal synteny scores are then averaged across all species used as reference species for OG synteny.
        SYN = statistics.mean(max_synteny_for_all_refSpecies)
        orthoGroupsMetrics[og]['SYN'] = SYN
    else:
        # Setting SYNTENY to NA for all OGs that were not represented by Reference Species genes.
        orthoGroupsMetrics[og]['SYN'] = 'NA'

    # return orthoGroupsMetrics
