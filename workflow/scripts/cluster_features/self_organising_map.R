library(kohonen)

set.seed(12345)

THREADS <- snakemake@config[['MAX_THREADS']]
n_x = 10
n_y = 10
n_hclust = 10

merged_orthogroup_features <- read.delim(snakemake@input[['features']])
df <- na.omit(merged_orthogroup_features[,-1])
sdf <- scale(df, center=T, scale=T)

pc <- princomp(sdf)
data_train_matrix <- as.matrix(pc$scores[,1:10])

som_grid <- somgrid(xdim = n_x, ydim = n_y, topo="hexagonal", toroidal=FALSE)
som_model <- supersom(data_train_matrix,
                 grid=som_grid,
                 rlen=300,
                 alpha=c(0.05,0.01),
                 keep.data = TRUE, cores=THREADS)

## SOM plot

pdf(file=snakemake@output[['plot']])

plot(som_model, type="changes") # Plotting training process convergence
plot(table(som_model$unit.classif), xlab='SOM cell', ylab='mapped orthogroups')
plot(som_model, type="count", shape="straight")
plot(som_model, type="dist.neighbours", shape="straight", main = "SOM neighbour distances")
plot(som_model, type="codes", shape="straight", palette.name = rainbow, main = "SOM components")

## use hierarchical clustering to cluster the codebook vectors
plot(som_model, type="codes", shape="straight", palette.name = rainbow, main = "SOM components")
som.hc <- cutree(hclust(object.distances(som_model, "codes"), method='ward.D2'), n_hclust)
add.cluster.boundaries(som_model, som.hc, lwd=4, col='blue')

dev.off()


## SOM cell memberships

orthogroups <- na.omit(merged_orthogroup_features)[,1]
cluster_ids <- som_model$unit.classif
memberships <- cbind(orthogroups, cluster_ids)

write.table(memberships, snakemake@output[['som_clusters']], quote = F, row.names = F, col.names = F, sep='\t')


## SOM supercell memberships

superclusters_ids <- c()
for(i in 1:nrow(memberships)){
  cluster_id <- as.integer(memberships[i, 'cluster_ids'])
  superclusters_ids <- c(superclusters_ids, som.hc[cluster_id])
}

supermemberships <- cbind(orthogroups, superclusters_ids)

write.table(supermemberships, snakemake@output[['som_hc_superclusters']], quote = F, row.names = F, col.names = F, sep='\t')
