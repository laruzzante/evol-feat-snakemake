library(kohonen)

THREADS <- snakemake@config[['MAX_THREADS']]

merged_orthogroup_features <- read.delim(snakemake@input[['features']])

data_train <- na.omit(merged_orthogroup_features[,-1])


ndims <- length(colnames(data_train))

data_train_matrix <- as.matrix(scale(data_train, center=TRUE, scale=TRUE))

som_grid <- somgrid(xdim = 25, ydim = 25, topo="hexagonal", toroidal=TRUE)

som_model <- supersom(data_train_matrix,
                 grid=som_grid,
                 rlen=1000,
                 alpha=c(0.05,0.01),
                 keep.data = TRUE, cores=THREADS)


pdf(file=snakemake@output[['plot']])

plot(som_model, type="changes")
plot(som_model, type="count", shape="straight")
plot(som_model, type="dist.neighbours", shape="straight", main = "SOM neighbour distances")
plot(som_model, type="codes", shape="straight", palette.name = rainbow, main = "SOM components")

## use hierarchical clustering to cluster the codebook vectors
plot(som_model, type="codes", shape="straight", palette.name = rainbow, main = "SOM components")
som.hc <- cutree(hclust(object.distances(som_model, "codes")), 10)
add.cluster.boundaries(som_model, som.hc, lwd=3, col='darkred')

dev.off()
