#!/usr/bin/env Rscript
library(topGO, quietly = T)

setwd('~/evol-feat_test/moz2_GO/')

go <- read.table('AGAP2GOS.txt', header = F)
geneID2GO <- readMappings(file = "AGAP2GOSf-topGO.txt")
geneUniverse <- names(geneID2GO)
genesOfInterest = read.table('~/evol-feat/dash/genesOfInterest.tsv', header=F)
genesOfInterest <- as.character(genesOfInterest$V1)
print(genesOfInterest)

geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse

myGOdata <- new("topGOdata", description="Evol-Feat SOM cell", ontology="BP", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
sg <- sigGenes(myGOdata)
str(sg)
numSigGenes(myGOdata)

resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
resultFisher
allRes <- GenTable(myGOdata, classicFisher=resultFisher, orderBy="resultFisher", ranksOf="classicFisher", topNodes=10)
allRes

length(usedGO(myGOdata))
