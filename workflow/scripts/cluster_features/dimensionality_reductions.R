library(uwot)
library(Rtsne)

set.seed(1234)

THREADS <- snakemake@config[['MAX_THREADS']]

merged_orthogroup_features <- read.delim(snakemake@input[[1]])
df <- na.omit(merged_orthogroup_features[,2:ncol(merged_orthogroup_features)])

## Scaling and centering
sdf <- scale(df, center = TRUE, scale = TRUE)

## PCA
pc <- princomp(sdf)
pc.coords <- pc$scores

## tSNE
perpl <- sqrt(nrow(pc.coords))
tsne <- Rtsne(pc.coords, perplexity = perpl, check_duplicates = FALSE, pca = FALSE, dims = 2, num_threads = THREADS)
tsne.coords <- tsne$Y

## UMAP
umap.coords <- umap(pc.coords, n_components = 2, n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = THREADS, metric = 'correlation')

## Writing output coordinate tables
write.table(pc.coords, snakemake@output[[1]], quote = F, row.names = F, col.names = F, sep='\t')
write.table(tsne.coords, snakemake@output[[2]], quote = F, row.names = F, col.names = F, sep='\t')
write.table(umap.coords, snakemake@output[[3]], quote = F, row.names = F, col.names = F, sep='\t')
