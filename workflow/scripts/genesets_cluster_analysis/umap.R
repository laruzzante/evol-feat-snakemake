input_features <- na.omit(read.delim(snakemake@input[[1]]))

library(uwot)

n_principal_components <- 10
pc_mean_sets <- prcomp(input_features[,-c(1,2)], center=TRUE, scale.=TRUE)$x[,1:n_principal_components]
train <- cbind(input_features[,c(1,2)], pc_mean_sets)

set.seed(1234)
tsne_on_pc <- Rtsne(train[,-c(1,2)], dims = 2, perplexity=30, verbose=TRUE, max_iter=500, check_duplicates=FALSE)
tsne <- Rtsne(input_features[,-c(1,2)], dims = 2, perplexity=30, verbose=TRUE, max_iter=500, check_duplicates=FALSE)


umap_on_pc <- umap(train[,-c(1,2), n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8)
umap <- umap(train[,-c(1,2), n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8, metric = 'correlation')

plot(umapdf, pch=19, col=t.blue, cex=0.2)
