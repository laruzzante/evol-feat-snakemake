library(topGO, quietly = T)
library(Rgraphviz, quietly = T)

go_universe <- snakemake@input[['go_universe']]
som_clusters <- snakemake@input[['som_clusters']]
ontology <- snakemake@params[['ontology']]

print(ontology)

output_table <- snakemake@output[['som_clusters_go']]
output_pdf <- snakemake@output[['som_clusters_go_dag']]

orthogroups2go <- readMappings(file = go_universe)
orthogroupsUniverse <- names(orthogroups2go)

str(head(orthogroups2go))

orthogroups_of_interest <- read.table(som_clusters)

check_and_delete <- function(appended_output_file){
  if (file.exists(appended_output_file)) {
    #Delete file if it exists
    file.remove(appended_output_file)
  }
}

check_and_delete(output_table)
check_and_delete(output_pdf)


pdf(output_pdf)

for(cluster in unique(sort(orthogroups_of_interest$V2))){

  myInterestingOGs <- as.character(orthogroups_of_interest[orthogroups_of_interest$V2==cluster,'V1'])

  print(paste0('Cluster ', cluster))
  print(myInterestingOGs)

  ogList <- factor(as.integer(orthogroupsUniverse %in% myInterestingOGs))
  names(ogList) <- orthogroupsUniverse
  str(ogList)

  GOdata <- new("topGOdata", ontology = ontology, allGenes = ogList,
                annot = annFUN.gene2GO, gene2GO = orthogroups2go)


  # In this scenario, when only a list of interesting genes is
  # provided, the user can use only tests statistics that are based on gene counts, like Fisher's exact test, Z score
  # and alike.


  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  print(resultFisher)

  resultFisher.weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  print(resultFisher.weight01)

  allRes <- GenTable(GOdata, classicFisher = resultFisher, weight01Fisher = resultFisher.weight01,
                     orderBy =  "weight01Fisher", ranksOf = "classicFisher", topNodes = 30)
  print(allRes)


  cat(paste0("Cluster ",cluster,"\n"), file = output_table, append = T)
  write.table(allRes, file=output_table, quote = F, sep='\t', col.names = T, row.names = F, append=T)

  showSigOfNodes(GOdata, score(resultFisher.weight01), firstSigNodes = 10, useInfo = 'all')

}

dev.off()
