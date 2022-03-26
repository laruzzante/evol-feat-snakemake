library(kohonen)
library(colorRamps)
library(GO.db)

setwd('~/evol-feat/')
data <- read.table(file='output_data/Insecta_OGs_AGAMPfiltered_Genes_EvolMetrics.tsv', header = T)


## BUILDING SOM
data_train <- data[,-c(1,2,3)]
rownames(data_train) <- data$gene_id

data_train_matrix <- as.matrix(scale(data_train))
som_grid <- somgrid(xdim=25, ydim=25, topo='hexagonal',
                    neighbourhood.fct='gaussian', toroidal=T)
sommap <- som(data_train_matrix,
              grid=som_grid,
              rlen=1000,
              alpha=c(0.05, 0.01),
              keep.data = TRUE)


## MAPPING SOM CELLS to GENE IDs AND OGs
gene_id <- row.names(data_train)
som_cell <- c()
for(i in 1:length(sommap$unit.classif)){
  som_cell <- c(som_cell, sommap$unit.classif[i])
}
gene_ids_2_sommap <- cbind(gene_id, som_cell)

gene_2_og <- data[,c('gene_id', 'pub_og_id')]
gene_2_og_2_sommap <- merge(gene_2_og, gene_ids_2_sommap, by = 'gene_id')


## MAPPING SOM CELLS, GENE IDs AND OGs to GO TERMS
AGAP2GOS <- read.delim("moz2_GO/AGAP2GOS.txt", header=FALSE)
colnames(AGAP2GOS) <- c('GO.term', 'Var2', 'gene_id')
df <- merge(gene_2_og_2_sommap, AGAP2GOS, by = 'gene_id')
df$som_cell <- as.numeric(df$som_cell)


## MAPPING SOM CELLS, GENE IDs, OGs AND GO TERMS to GO TERM ANNOTATIONS
GO_ids <- df$GO.term
GO_annotations <- c()
for(i in 1:length(GO_ids)){
  GO_annotations <- c(GO_annotations, as.character(Term(as.character(GO_ids[i]))))
}
df$GO.annotations <- GO_annotations


## WRITING OUTPUT FILE
write.table(df, file = 'moz2_GO/AGAP2GOS_2_somcells.tsv',
            row.names = F, sep = '\t', quote = F)

