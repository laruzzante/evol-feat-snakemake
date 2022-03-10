

## Scaling and centering
sdf <- scale(df, center = TRUE, scale = TRUE)

## PCA
pc <- princomp(na.omit(sdf))

## tSNE
perpl <- sqrt(nrow(sdf))
tsne <- Rtsne(sdf, perplexity = perpl, check_duplicates=FALSE, pca = FALSE, dims = 2, num_threads = 8)

## UMAP
umap <- umap(sdf, n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8, metric = 'correlation')
