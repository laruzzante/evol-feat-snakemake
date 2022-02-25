for og in sorted(ogDict.keys()):
    for spec in all_species:
        ogCountsDict[og][spec] = 0

    for specInOG, genes in sorted(ogDict[og].items()):
        ogCountsDict[og][specInOG] += len(genes)

for og, species in sorted(ogCountsDict.items()):
    outline = str(len(ogDict[og].keys())) + '\t' + og
    for spec in sorted(species):
        outline = outline + '\t' + str(ogCountsDict[og][spec])
    print(outline)
