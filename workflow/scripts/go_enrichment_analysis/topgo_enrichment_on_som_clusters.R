library(topGO)

go_universe <- snakemake@input[['go_universe']]
som_clusters <- snakemake@input[['som_clusters']]

orthogroups2go <- readMappings(file = go_universe)
orthogroupsUniverse <- names(orthogroups2go)

orthogroups_of_interest <- read.table(som_clusters)

som_clusters_go <- c()
for cluster in unique(orthogroups_of_interest$V2){
  orthogroups_per_cluster <- as.character(orthogroups_of_interest[orthogroups_of_interest$V2==cluster,'V1'])
  print(paste0('Cluster ', cluster))
  print(orthogroups_per_cluster)

  ogList <- factor(as.integer(orthogroupsUniverse %in% orthogroups_per_cluster))
  names(ogList) <- orthogroupsUniverse

  myGOdata <- new(paste0("topGOdata",clsuter), description="Evol-Feat SOM cell", ontology="BP", allGenes=ogList, annot=annFUN.gene2GO, gene2GO=orthogroups2go)
  sg <- sigGenes(myGOdata)
  str(sg)
  numSigGenes(myGOdata)

  resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
  print(resultFisher)
  allRes <- GenTable(myGOdata, classicFisher=resultFisher, orderBy="resultFisher", ranksOf="classicFisher", topNodes=10)
  print(allRes)
  som_clusters_go <- c(som_clusters_go, allRes)

  print(length(usedGO(myGOdata)))
}

write.table(som_clusters_go, snakemake@output[['som_clusters_go']], quote = F, row.names = F, col.names = F, sep='\t')
