library(kohonen)

set.seed(12345)

THREADS <- snakemake@config[['MAX_THREADS']]
n_x = 10
n_y = 10
n_supercl = 12 # For distinct colours visualisation purpose, we will not compute more than 12 superclusters

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
plot(som_model, type="mapping", shape="straight", palette.name = rainbow, main = "SOM components")


## use hierarchical clustering to cluster the codebook vectors
# som.hc <- cutree(hclust(pearson(object.distances(som_model, "codes")), method='average'), n_supercl)
# som.hc <- cutree(hclust(object.distances(som_model, "codes"), method='ward.D2'), n_supercl)
# som.kmeans <- kmeans(object.distances(som_model, "codes"), centers = n_supercl)

# library(dbscan)
# som.dbscan <- dbscan(object.distances(som_model, "codes"), eps = 0.15, minPts = 2)
# som.hdbscan <- hdbscan(object.distances(som_model, "codes"), minPts = 2)

codes <- do.call(rbind.data.frame, som_model$codes)
som.kmeans <- kmeans(codes, centers = n_supercl)

## COMPUTING


library(RColorBrewer)
# display.brewer.all() ## Decomment to show all available palettes.
supercl <- som.kmeans$cluster

colours <- brewer.pal(n_supercl, 'Paired') ## Sticking to Paired palette as its the one with the most distintive colours (12)
colours <- cbind(as.character(1:n_supercl), colours)
colnames(colours) <- c('set', 'colour')
colours <- data.frame(colours)
supercl <- data.frame(set = as.character(supercl))

library(dplyr)
som.supercl.colours <- inner_join(supercl, colours)

plot(som_model, type="codes", shape="straight", palette.name = rainbow, main = "SOM components", bgcol=som.supercl.colours$colour)
add.cluster.boundaries(som_model, supercl$set, lwd=2, col='black')
legend(x=-1.5, y = 8, legend = colours$set, fill = colours$colour, title = 'super_cl')


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
  superclusters_ids <- c(superclusters_ids, clusters[cluster_id])
}

supermemberships <- cbind(orthogroups, superclusters_ids)

write.table(supermemberships, snakemake@output[['som_superclusters']], quote = F, row.names = F, col.names = F, sep='\t')
