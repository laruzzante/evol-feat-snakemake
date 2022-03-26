source("~/agambiae-immunity-families/livio/scripts/functions.R")
library(rgl)
options(rgl.printRglwidget = TRUE)
library(Gmedian)
library(dbscan)
library(Rtsne)
library(colorRamps)
library(dendextend)
library(uwot)

library(factoextra)



# Transparent blue
t.blue <- rgb(0, 0, 255, max = 255, alpha = 25, names = "blue50")

# Succession of N most distinct colours
n <- 300
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colours = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


setwd("~/evol-feat-snakemake/workflow/output")

merged_orthogroup_features <- read.delim("~/evol-feat-snakemake/workflow/output/merged_orthogroup_features.tsv")

df <- na.omit(merged_orthogroup_features[,2:ncol(merged_orthogroup_features)])

sdf <- scale(df)


## PCA
pc <- princomp(sdf)

plot(pc$scores[,1],pc$scores[,2], pch=19, col=t.blue, cex=0.2)

print(pc$loadings)

summary(pc)
# Control variable colors using their contributions
fviz_pca_var(pc, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) +  theme_minimal()

fviz_pca_biplot(pc, pointsize=0.2, label='var', repel=T, col.var = 'red', col.ind='blue', alpha.ind=0.2) + theme_minimal()
plot3d(pc$scores[,3], y=pc$scores[,1], z=pc$scores[,2],col=t.blue)

## tSNE

max_perplexity <- floor((nrow(sdf) - 1) / 3)
perpl <- sqrt(nrow(sdf))
tsnedf <- Rtsne(pc$scores, perplexity = perpl, check_duplicates=FALSE, pca = FALSE, dims = 2, num_threads = 6) # Defaults: perplexity=30 and pca=TRUE ... i.e. we have to specify that we do not want a PCA beforehand
# plot3d(tsnedf$Y,col=t.blue) # Only for dims=3
plot(tsnedf$Y, pch=19, col=t.blue, cex=0.2)

n <- 10

## UMAP
# umapdf <- umap(sdf, n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8, metric = 'correlation')
umapdf <- umap(pc$scores, n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8, metric = 'correlation') # Default settings for everything apart from min_dist, set from 0.01 to 0.001 and metric from Euclidean to correlation
# min_dist
# The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on
# the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set relative to the
# spread value, which determines the scale at which embedded points will be spread out.

plot(umapdf, pch=19, col=t.blue, cex=0.2)

## KMEANS

cl.kmeans <- kmeans(sdf, n)
plot(pc$scores[,1],pc$scores[,2], pch=19, col=cl.kmeans$cluster, cex=0.2)
plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], pch=19, col=cl.kmeans$cluster, cex=0.2)

cl.kmeans <- kmeans(pc$scores[,1:3], n)
plot(pc$scores[,1],pc$scores[,2], pch=19, col=cl.kmeans$cluster, cex=0.2)
plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], pch=19, col=cl.kmeans$cluster, cex=0.2)

cl.kmeans <- kmeans(pc$scores, n)
plot(pc$scores[,1],pc$scores[,2], pch=19, col=cl.kmeans$cluster, cex=0.2)
plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], pch=19, col=cl.kmeans$cluster, cex=0.2)

## KMEDIANS

cl.kmedian <- kGmedian(sdf, ncenters=n, iter.max = 100)
plot(pc$scores[,1],pc$scores[,2], pch=19, col=cl.kmedian$cluster, cex=0.2)
plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], pch=19, col=cl.kmedian$cluster, cex=0.2)

cl.kmedian <- kGmedian(pc$scores, ncenters=n, iter.max = 100)
plot(pc$scores[,1],pc$scores[,2], pch=19, col=cl.kmedian$cluster, cex=0.2)
plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], pch=19, col=cl.kmedian$cluster, cex=0.2)

### OPTICS

## PCA OPTICS
cl.optics <- optics(pc$scores)
cl.optics.cut <- extractDBSCAN(cl.optics, eps_cl = 1)
plot(pc$scores[,1],pc$scores[,2], pch=1, col=cl.optics.cut$cluster+1L, cex=0.2)
plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], radius=0.2, col=cl.optics.cut$cluster+1L)

dbscan::kNNdistplot(tsnedf$Y, k = 5)
abline(h = 0.3, lty = 2) # The knee seems to be around 0.15, hence we plot a line at h=0.15 just to see the actual intersection
cl.optics <- optics(tsnedf$Y, eps = 0.3, minPts = 5)
cl.optics.cut <- extractDBSCAN(cl.optics, eps_cl = 0.3)
palette <- c()
for(cl in cl.optics.cut$cluster){
  x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
  while(x > length(colours)){
    x <- x - length(colours)
  }
  palette <- c(palette, colours[x])
}
plot(tsnedf$Y, pch=1, col=palette, cex=0.2)
plot3d(tsnedf$Y, radius=0.2, col=palette)


cl.optics <- optics(umapdf, eps = 0.15, minPts = 10)
cl.optics.cut <- extractDBSCAN(cl.optics, eps_cl = 0.15)
palette <- c()
symbols <- c()
for(cl in cl.optics.cut$cluster){
  if(cl == 0){
    palette <- c(palette, '#000000') # When cluster == 0, i.e. points that could not be assigned to any cluster, then colour is black
    symbols <- c(symbols, 4)# When cluster == 0, i.e. points that could not be assigned to any cluster, then symbol is a 'x'
  } else {
    x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
    while(x > length(colours)){
      x <- x - length(colours)
    }
    palette <- c(palette, colours[x])

    symbol <- cl%%25
    if(symbol == 4){
      symbol <- 0
    }
    symbols <- c(symbols, symbol)
  }
}
plot(umapdf, pch=symbols, col=palette, cex=0.2, lwd=0.2)
plot3d(umapdf, radius=0.2, col=palette)

### DBSCAN

## DBSCAN on PCA
dbscan::kNNdistplot(pc$scores, k = 10) ## By plotting this one we can check where the knee is. Careful, k must be equal to minPoints used above in dbscan
abline(h = 1.5, lty = 2) # The knee seems to be around 0.15, hence we plot a line at h=0.15 just to see the actual intersection

cl.dbscan <- dbscan(pc$scores, eps=1.5, weights = pc$sdev.^2 / sum(pc$sdev), minPts = 10) # I have added weights on the clustering, i.e. each PC is being weighted by the
n_clusters <- length(unique(cl.dbscan$cluster))
palette <- c()
symbols <- c()
for(cl in cl.dbscan$cluster){
  if(cl == 0){
    palette <- c(palette, '#000000') # When cluster == 0, i.e. points that could not be assigned to any cluster, then colour is black
    symbols <- c(symbols, 4)# When cluster == 0, i.e. points that could not be assigned to any cluster, then symbol is a 'x'
  } else {
    x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
    while(x > length(colours)){
      x <- x - length(colours)
    }
    palette <- c(palette, colours[x])

    symbol <- cl%%25
    if(symbol == 4){
      symbol <- 0
    }
    symbols <- c(symbols, symbol)
  }
}
mapping <- list('palette'=palette, 'symbols'=symbols)

plot(pc$scores[,c(1,2)], pch=unlist(mapping['symbols']), col=unlist(mapping['palette']), cex=0.2, lwd=0.2)

orthogroups <- na.omit(merged_orthogroup_features)[,1]
memberships <- cl.dbscan$cluster

memerbships <- cbind(orthogroups, memberships)

# coords <- Rtsne(pc$scores[,1:5], perplexity = perpl, check_duplicates=FALSE, pca = FALSE, dims = 2, num_threads = 6) # Defaults: perplexity=30 and pca=TRUE ... i.e. we have to specify that we do not want a PCA beforehand
# plot(coords$Y, pch=1, col=cl.dbscan$cluster+1L, cex=0.2)
# plot(pc$scores, pch=1, col=cl.dbscan$cluster+1L, cex=0.2)
plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], radius=3, col=cl.dbscan$cluster+1L)

# DBSCAN on tSNE (tSNE might have been done on PCA values, depending on how it was called before)
cl.dbscan <- dbscan(tsnedf$Y, eps=0.3, minPts = 2)
cl.dbscan_filtered <- subset(cl.dbscan, cl.dbscan$cluster != 0)
n_clusters <- length(unique(cl.dbscan$cluster))
palette <- c()
for(cl in cl.dbscan$cluster){
  x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
  while(x > length(colours)){
    x <- x - length(colours)
  }
  palette <- c(palette, colours[x])
}
plot(tsnedf$Y, pch=1, col=palette, cex=0.2)
# plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], radius=3, col=cl.dbscan$cluster+1L)


# DBSCAN on UMAP

# Method for determining the optimal eps value
# The method proposed here consists of computing the he k-nearest neighbor distances in a matrix of points.
# The idea is to calculate, the average of the distances of every point to its k nearest neighbors. The value of k will be specified by the user and corresponds to MinPts.
# Next, these k-distances are plotted in an ascending order. The aim is to determine the “knee”, which corresponds to the optimal eps parameter.
# A knee corresponds to a threshold where a sharp change occurs along the k-distance curve.
# The function kNNdistplot() [in dbscan package] can be used to draw the k-distance plot:
dbscan::kNNdistplot(umapdf, k = 10) ## By plotting this one we can check where the knee is. Careful, k must be equal to minPoints used above in dbscan
abline(h = 0.15, lty = 2) # The knee seems to be around 0.15, hence we plot a line at h=0.15 just to see the actual intersection

cl.dbscan <- dbscan(umapdf, eps=0.15, minPts = 10)
n_clusters <- length(unique(cl.dbscan$cluster))
palette <- c()
symbols <- c()
for(cl in cl.dbscan$cluster){
  if(cl == 0){
    palette <- c(palette, '#000000') # When cluster == 0, i.e. points that could not be assigned to any cluster, then colour is black
    symbols <- c(symbols, 4)# When cluster == 0, i.e. points that could not be assigned to any cluster, then symbol is a 'x'
  } else {
    x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
    while(x > length(colours)){
      x <- x - length(colours)
    }
    palette <- c(palette, colours[x])

    symbol <- cl%%25
    if(symbol == 4){
      symbol <- 0
    }
    symbols <- c(symbols, symbol)
  }
}
plot(umapdf, pch=symbols, col=palette, cex=0.2, lwd=0.2)

library(ggplot2)
library(cowplot)
# Basic scatter plot
print(length(cl.hdbscan$cluster[cl.hdbscan$cluster==0])) ## Printing N of unassigned orthogroups (noise points)
print(length(unique(cl.hdbscan$cluster)) -1 ) ## Printing N of clusters, minus 1, i.e. to remove cluster id 0 (unassigned)
umap.df <-data.frame(umapdf)
ggplot(umap.df, aes(x=X1, y=X2)) +
  geom_point(shape=symbols, color=palette, fill=palette, size=0.8, alpha=0.2, stroke=0.4) +
  theme_bw() + xlab(label='UMAP.X') + ylab(label='UMAP.Y')

## HDBSCAN
HDBminPoints = 10
sset <- umapdf[sample(nrow(umapdf), 20000), ]
cl.hdbscan <- hdbscan(sset, minPts = HDBminPoints, gen_hdbscan_tree = T)
# hdb_filtered <- subset(hdb, hdb$cluster != 0)

palette <- c()
symbols <- c()
for(cl in cl.hdbscan$cluster){
  if(cl == 0){
    palette <- c(palette, '#000000') # When cluster == 0, i.e. points that could not be assigned to any cluster, then colour is black
    symbols <- c(symbols, 4)# When cluster == 0, i.e. points that could not be assigned to any cluster, then symbol is a 'x'
  } else {
    x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
    while(x > length(colours)){
      x <- x - length(colours)
    }
    palette <- c(palette, colours[x])

    symbol <- cl%%25
    if(symbol == 4){
      symbol <- 0
    }
    symbols <- c(symbols, symbol)
  }
}
plot(sset, pch=symbols, col=palette, cex=0.2, lwd=0.2)

library(ggplot2)
library(cowplot)
# Basic scatter plot
print(length(cl.dbscan$cluster[cl.dbscan$cluster==0])) ## Printing N of unassigned orthogroups (noise points)
print(length(unique(cl.dbscan$cluster)) -1 ) ## Printing N of clusters, minus 1, i.e. to remove cluster id 0 (unassigned)
umap.df <-data.frame(umapdf)
ggplot(umap.df, aes(x=X1, y=X2)) +
  geom_point(shape=symbols, color=palette, fill=palette, size=0.8, alpha=0.2, stroke=0.4) +
  theme_bw() + xlab(label='UMAP.X') + ylab(label='UMAP.Y')

plot(sset, col=hdb$cluster,
     pch=ifelse(hdb$cluster == 0, 8, 1), # Mark noise as star
     cex=ifelse(hdb$cluster == 0, 0.5, 0.75), # Decrease size of noise
     xlab=NA, ylab=NA)
colors <- sapply(1:length(hdb$cluster),
                 function(i) adjustcolor(palette()[(hdb$cluster+1)[i]], alpha.f = hdb$membership_prob[i]))
points(DS3, col=colors, pch=20)


## HIERARCHICAL CLUSTERING
bootstraps_families <- pvclust(t(pc$scores[1:5000,]), method.dist = 'euclidean', method.hclust = 'ward.D2', nboot = 100, parallel = T, iseed = 1234)
bootstraps_metrics <- pvclust(pc$scores[1:5000,], method.dist = 'euclidean', method.hclust = 'ward.D2', nboot = 100, parallel = T, iseed = 1234)
# hm <- Heatmap(pc$scores[1:1000,], cluster_rows = bootstraps_families$hclust, cluster_columns = bootstraps_metrics$hclust, heatmap_legend_param = list(title = "Scaled\nMetrics"))
# print(hm)

metrics <- as.dendrogram(bootstraps_metrics$hclust)
families <- as.dendrogram(bootstraps_families$hclust)

# AU p-values from pvclust
metrics_AU <- bootstraps_metrics$edges$au
families_AU <- bootstraps_families$edges$au


# Ordering the dendrogram nodes from pvclust following as.dendrogram node order, WHICH IS DIFFERENT!@#$%&*!!
# And as.dendrogram also includes leaves indices within nodes. So be careful that leaves and nodes are all mixed togheter.
# I recognize leaves by their height of 0.
# I'm using height of nodes as the mapping variable from one pvclust dendrogram to as.dendrogram becuase it's the only
# common factor I could use. So, pvclust dendrogram nodes are named by their height score and then ordered following the
# height vector from as.dendrogram which corresponds to the actual plotting order of the nodes. The limitation is that
# we must be careful in the rare eventuality that two nodes have the exact same height.
names(metrics_AU) <- bootstraps_metrics$hclust$height # metrics AU p-values
dendrogram_metrics_heights <- get_nodes_attr(metrics, "height", include_leaves = T)
metrics_AU_sorted <- c()
for(node_height in dendrogram_metrics_heights){
  if(node_height != 0.0) metrics_AU_sorted <- c(metrics_AU_sorted, metrics_AU[as.character(node_height)])
  else metrics_AU_sorted <- c(metrics_AU_sorted, NA)
}

names(families_AU) <- bootstraps_families$hclust$height # same for families
dendrogram_families_heights <- get_nodes_attr(families, "height", include_leaves = T)
families_AU_sorted <- c()
for(node_height in dendrogram_families_heights){
  if(node_height != 0.0) families_AU_sorted <- c(families_AU_sorted, families_AU[as.character(node_height)])
  else families_AU_sorted <- c(families_AU_sorted, NA)
}


## Heatmap

# Setting AU p-values colour-scale in dendrograms
AU_colorpalette <- c(colorRampPalette(brewer.pal(9, "Greens")[3:9])(9))
# image(1:length(AU_colorpalette), 1, as.matrix(1:length(AU_colorpalette)), col=AU_colorpalette)

map2color <- function(AU_list){
  AU_colors <- c()
  for(AU in AU_list){
    if(is.na(AU)) AU_colors <- c(AU_colors, NA)
    else if(AU < 0.6) AU_colors <- c(AU_colors, AU_colorpalette[1])
    else if(AU < 0.65) AU_colors <- c(AU_colors, AU_colorpalette[2])
    else if(AU < 0.7) AU_colors <- c(AU_colors, AU_colorpalette[3])
    else if(AU < 0.75) AU_colors <- c(AU_colors, AU_colorpalette[4])
    else if(AU < 0.8) AU_colors <- c(AU_colors, AU_colorpalette[5])
    else if(AU < 0.85) AU_colors <- c(AU_colors, AU_colorpalette[6])
    else if(AU < 0.9) AU_colors <- c(AU_colors, AU_colorpalette[7])
    else if(AU < 0.95) AU_colors <- c(AU_colors, AU_colorpalette[8])
    else if(AU <= 1.0) AU_colors <- c(AU_colors, AU_colorpalette[9])
    else print('Warning: AU value not in range ]-Inf,1]')
  }
  return(AU_colors)
}

metrics_colors_AU <- map2color(metrics_AU_sorted)
families_colors_AU <- map2color(families_AU_sorted)

heatmap_colors <- rev(unlist(strsplit('1E009F-2531C7-0056FF-3C8DF9-6BAAF0-BCDBF7-F8F987-FFCA00-F9A500-FF6500-FF3800-D20000', '-')))
heatmap_colors <- paste0("#", heatmap_colors)


create_bins <- function(hm){
  legend_bins <- c()
  for(i in seq(1, length(hm$colorTable$color))){
    if(i == 1){
      bin <- paste0('[', round(hm$colorTable$low[i], 2), ', ', round(hm$colorTable$high[i], 2), '[')
    } else if(i == length(hm$colorTable$color)) {
      bin <- paste0('[', round(hm$colorTable$low[i], 2), ', ', round(hm$colorTable$high[i], 2), ']')
    } else {
      bin <- paste0('[', round(hm$colorTable$low[i], 2), ', ', round(hm$colorTable$high[i], 2), '[')
    }
    legend_bins <- c(legend_bins, bin)
  }
  return(legend_bins)
}

Rowv <- families %>%
  #set("branches_k_color", k = 10) %>% set("branches_lwd", 3) %>%  # setting the number of separate colour-coded clades
  set("nodes_pch", 20) %>%  # node point type
  set("nodes_cex", 1.5) %>%  # node point size
  set("nodes_col", families_colors_AU, warn=T) %>% # node point color
  set("clear_leaves") #%>%

Colv <- metrics %>%
  #set("branches_k_color", k = 5) %>% set("branches_lwd", 3) %>%  # setting the number of separate colour-coded clades
  set("nodes_pch", 20) %>%  # node point type
  set("nodes_cex", 1.5) %>%  # node point size
  set("nodes_col", metrics_colors_AU, warn=T) %>% # node point color
  set("clear_leaves") #%>%

hm <- heatmap.2(pc$scores[1:5000,], Rowv = Rowv, Colv = Colv, trace = "none", labRow = FALSE,
                density.info = "none", col = rev(heatmap_colors), scale = "row", key = FALSE)
legend(x=0.01,y=0.5, title = '% AU Support',
       legend = rev(c("<60","60-65", "65-70", "70-75", "75-80", "80-85", "85-90", "90-95", "95-100")),
       col = rev(AU_colorpalette),
       lty= 1,
       lwd = 5,
       cex=0.5
)
legend <- legend(x=0,y=1, title = 'Scaled Metrics PCs',
       legend = rev(create_bins(hm)),
       col = heatmap_colors,
       lty= 1,
       lwd = 5,
       cex=0.5
)
