
library(pvclust)
library(ComplexHeatmap)

scale_features <- function(df){
  scaled_df <- df
  for(i in 1:ncol(df)){
    if(sapply(df[i], is.numeric)){
      scaled_df[i] <- scale(df[i])
    }
  }
  return(scaled_df)
}

format_features_matrix <- function(df, aggregator_function=aggregator_function){
  df_aggregate <- aggregate(df[,3:ncol(df)], by=list(df[['set']]), FUN = aggregator_function, na.rm=TRUE)
  rownames(df_aggregate) <- df_aggregate$Group.1
  evol_features_matrix <- as.matrix(df_aggregate[2:ncol(df_aggregate)])
  return(evol_features_matrix)
}

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

plot_dendrogram_heatmap <- function(df, nboot=10000, aggregator_function, method.dist=pearson, method.hclust="average", iseed=1234){
  scaled_df <- scale_features(df)
  evol_features_matrix <- format_features_matrix(scaled_df, aggregator_function)

  bootstraps_features <- pvclust(evol_features_matrix, method.dist = method.dist, method.hclust = method.hclust, nboot = nboot, parallel = T, iseed = iseed)
  plot(bootstraps_features)
  pvrect(bootstraps_features)

  bootstraps_sets <- pvclust(t(evol_features_matrix), method.dist = method.dist, method.hclust = method.hclust, nboot = nboot, parallel = T, iseed = iseed)
  plot(bootstraps_sets)
  pvrect(bootstraps_sets)

  # ComplexHeatmap package used to combine custom pvclust dendrograms with heatmap
  hm <- Heatmap(evol_features_matrix, cluster_rows = bootstraps_sets$hclust, cluster_columns = bootstraps_features$hclust, heatmap_legend_param = list(title = "Scaled\nFeatures"))
  print(hm)
}
