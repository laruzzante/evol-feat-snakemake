source("~/agambiae-immunity-families/livio/scripts/functions.R")
library(rgl)
options(rgl.printRglwidget = TRUE)
library(Gmedian)
library(dbscan)
library(Rtsne)
library(colorRamps)
library(dendextend)
library(uwot)

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
pc <- princomp(na.omit(sdf))

plot(pc$scores[,1],pc$scores[,2], pch=19, col=t.blue, cex=0.2)

print(pc$loadings)

summary(pc)

plot3d(pc$scores[,3], y=pc$scores[,1], z=pc$scores[,2],col=t.blue)

## tSNE

max_perplexity <- floor((nrow(sdf) - 1) / 3)
tsnedf <- Rtsne(sdf, perplexity = 30, check_duplicates=FALSE, pca = TRUE, dims = 2, num_threads = 6) # Defaults: perplexity=30 and pca=TRUE ... i.e. we have to specify that we do not want a PCA beforehand
# plot3d(tsnedf$Y,col=t.blue) # Only for dims=3
plot(tsnedf$Y, pch=19, col=t.blue, cex=0.2)

n <- 10

## UMAP
umapdf <- umap(sdf, n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8, metric = 'correlation')
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

## OPTICS

cl.optics <- optics(pc$scores)
cl.optics.cut <- extractDBSCAN(cl.optics, eps_cl = 1)
plot(pc$scores[,1],pc$scores[,2], pch=1, col=cl.optics.cut$cluster+1L, cex=0.2)
plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], radius=0.2, col=cl.optics.cut$cluster+1L)

cl.optics <- optics(tsnedf$Y, eps = 0.4)
cl.optics.cut <- extractDBSCAN(cl.optics, eps_cl = 0.4)
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


cl.optics <- optics(umapdf, eps = 0.3)
cl.optics.cut <- extractDBSCAN(cl.optics, eps_cl = 0.3)
palette <- c()
for(cl in cl.optics.cut$cluster){
  x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
  while(x > length(colours)){
    x <- x - length(colours)
  }
  palette <- c(palette, colours[x])
}
plot(umapdf, pch=1, col=palette, cex=0.2)
plot3d(umapdf, radius=0.2, col=palette)

## DBSCAN

# DBSCAN on PCA
cl.dbscan <- dbscan(pc$scores, eps=0.3, weights = pc$sdev.^2 / sum(pc$sdev), minPts = 2) # I have added weights on the clustering, i.e. each PC is being weighted by the
plot(tsnedf$Y, pch=1, col=cl.dbscan$cluster+1L, cex=0.2)
plot3d(x=pc$scores[,1],y=pc$scores[,2],z=pc$scores[,3], radius=3, col=cl.dbscan$cluster+1L)

# DBSCAN on tSNE (tSNE might have been done on PCA values, depending on how it was called before)
cl.dbscan <- dbscan(tsnedf$Y, eps=0.3, minPts = 2)
# amount of proportional variance it explains
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
cl.dbscan <- dbscan(umapdf, eps=0.3, minPts = 2) # I have added weights on the clustering, i.e. each PC is being weighted by the
# amount of proportional variance it explains
n_clusters <- length(unique(cl.dbscan$cluster))
palette <- c()
for(cl in cl.dbscan$cluster){
  x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
  while(x > length(colours)){
    x <- x - length(colours)
  }
  palette <- c(palette, colours[x])
}
plot(umapdf, pch=1, col=palette, cex=0.2)

# Method for determining the optimal eps value
# The method proposed here consists of computing the he k-nearest neighbor distances in a matrix of points.
# The idea is to calculate, the average of the distances of every point to its k nearest neighbors. The value of k will be specified by the user and corresponds to MinPts.
# Next, these k-distances are plotted in an ascending order. The aim is to determine the “knee”, which corresponds to the optimal eps parameter.
# A knee corresponds to a threshold where a sharp change occurs along the k-distance curve.
# The function kNNdistplot() [in dbscan package] can be used to draw the k-distance plot:
dbscan::kNNdistplot(pc$scores, k = 10) ## By plotting this one we can check where the knee is. Careful, k must be equal to minPoints used above in dbscan
abline(h = 2, lty = 2) # The knee seems to be around 2, hence we plot a line at h=2 just to see the actual intersection


## HDBSCAN

cl.hdbscan <- hdbscan(pc$scores, minPts = 25)
plot(pc$scores, col=cl.hdbscan$cluster, 
     pch=ifelse(cl.hdbscan$cluster == 0, 8, 1), # Mark noise as star
     cex=ifelse(cl.hdbscan$cluster == 0, 0.5, 0.75), # Decrease size of noise
     xlab=NA, ylab=NA)
colors <- sapply(1:length(cl2$cluster), 
                 function(i) adjustcolor(palette()[(cl2$cluster+1)[i]], alpha.f = cl2$membership_prob[i]))
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
