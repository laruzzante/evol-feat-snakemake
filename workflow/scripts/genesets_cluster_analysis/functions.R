library(pvclust, quietly=TRUE)
library(ComplexHeatmap, quietly=TRUE)

scale_features <- function(df){
  scaled_df <- df
  for(i in 1:ncol(df)){
    if(sapply(df[i], is.numeric)){
      scaled_df[i] <- scale(df[i], center=TRUE, scale=TRUE)
    }
  }
  return(scaled_df)
}

# Averaging gene features per gene set, need to specify function: mean, median, or other
format_features_matrix <- function(df, aggregator_function='mean'){
  df_aggregate <- aggregate(df[,3:ncol(df)], by=list(df[['set']]), FUN=aggregator_function, na.rm=TRUE)
  rownames(df_aggregate) <- df_aggregate$Group.1
  evol_features_matrix <- as.matrix(df_aggregate[2:ncol(df_aggregate)])
  return(evol_features_matrix)
}

# Defining correlation functions for pvclust dist calling method
pearson <- function(x){
  x <- as.matrix(x)
  res <- as.dist(1-cor(x, method="pearson"))
  attr(res, "method") <- "pearson"
  return(res)
}
spearman <-  function(x){
  x <- as.matrix(x)
  res <- as.dist(1-cor(x, method="spearman"))
  attr(res, "method") <- "spearman"
  return(res)
}
kendall <- function(x){
  x <- as.matrix(x)
  res <- as.dist(1-cor(x, method="kendall"))
  attr(res, "method") <- "kendell"
  return(res)
}

plot_pvclust_heatmap <- function(matrix, nboot=10000, aggregator_function='mean', method.dist='euclidean', method.hclust='ward.d2', iseed=1234){
  sets_dendrogram <- pvclust(t(matrix), method.dist=method.dist, method.hclust=method.hclust, nboot=nboot, parallel=T, iseed=iseed)
  features_dendrogram <- pvclust(matrix, method.dist=method.dist, method.hclust=method.hclust, nboot=nboot, parallel=T, iseed=iseed)
  # ComplexHeatmap package used to combine custom pvclust dendrograms with heatmap
  title <- paste('distance:', deparse(substitute(method.dist)),
                 '| method:', deparse(substitute(method.hclust)),
                 '\nsets summary:', aggregator_function, '\nboostraps:', nboot)
  heatmap <- Heatmap(matrix, cluster_rows=sets_dendrogram$hclust, cluster_columns=features_dendrogram$hclust,
                     column_title=title, heatmap_legend_param=list(title="Scaled\nFeatures"))

  print(heatmap)
  plot(sets_dendrogram)
  pvrect(sets_dendrogram)
  plot(features_dendrogram)
  pvrect(features_dendrogram)

  clusters <- c(heatmap, sets_dendrogram, features_dendrogram)
  names <- c('heatmap', 'sets_dendrogram', 'features_dendrogram')
  heatmap_clusters <- list('clusters' = clusters, 'names'=names)
  return(heatmap_clusters)
}
