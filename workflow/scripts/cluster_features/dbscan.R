# DBSCAN on Principal Components
cl.dbscan <- dbscan(pc$scores, eps=0.3, weights = pc$sdev.^2 / sum(pc$sdev), minPts = 2) # I have added weights on the clustering, i.e. each PC is being weighted by the
plot(tsnedf$Y, pch=1, col=cl.dbscan$cluster+1L, cex=0.2)


# DBSCAN on tSNE (tSNE might have been done on PCA values, depending on how it was called before)
cl.dbscan <- dbscan(tsnedf$Y, eps=0.3, minPts = 2) # I have added weights on the clustering, i.e. each PC is being weighted by the
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
