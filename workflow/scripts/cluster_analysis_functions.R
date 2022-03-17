library(pvclust, quietly=TRUE)
suppressPackageStartupMessages(library(ComplexHeatmap, quietly=TRUE))
library(RColorBrewer)


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


plot_pvclust_basic_heatmap <- function(matrix, nboot=100, scaling='columns', method.dist='euclidean', method.hclust='ward.d2', iseed=1234){
  orthogroups_dendrogram <- pvclust(t(matrix), method.dist=method.dist, method.hclust=method.hclust, nboot=nboot, parallel=T, iseed=iseed)
  features_dendrogram <- pvclust(matrix, method.dist=method.dist, method.hclust=method.hclust, nboot=nboot, parallel=T, iseed=iseed)
  # ComplexHeatmap package used to combine custom pvclust dendrograms with heatmap
  title <- paste('distance:', deparse(substitute(method.dist)),
                 '| method:', deparse(substitute(method.hclust)),
                 '\nscaling by:', scaling, '\nboostraps:', nboot)
  heatmap <- Heatmap(matrix, cluster_rows=orthogroups_dendrogram$hclust, cluster_columns=features_dendrogram$hclust,
                     column_title=title, show_row_names=FALSE, heatmap_legend_param=list(title="Scaled\nFeatures"))

  print(heatmap)
  plot(features_dendrogram)
  pvrect(features_dendrogram)

  clusters <- c(heatmap, orthogroups_dendrogram, features_dendrogram)
  names <- c('heatmap', 'sets_dendrogram', 'features_dendrogram')
  heatmap_clusters <- list('clusters' = clusters, 'names'=names)
  return(heatmap_clusters)
}


plot_pvclust_basic_heatmap_on_sets <- function(matrix, nboot=10000, aggregator_function='mean', method.dist='euclidean', method.hclust='ward.d2', iseed=1234){
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


map_palette_to_clusters <- function(cluster.object){

  ## Definine colour palette for plots
  # Succession of N most distinct colours
  n <- 300
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  colours = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


  palette <- c()
  symbols <- c()
  for(cl in cluster.object$cluster){
    if(cl == 0){
      palette <- c(palette, '#000000') # When cluster == 0, i.e. points that could not be assigned to any cluster, then colour is black
      symbols <- c(symbols, 4)# When cluster == 0, i.e. points that could not be assigned to any cluster, then symbol is a 'x'
    } else {
      x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
      while(x > length(colours)){
        x <- x - length(colours)
      }
      palette <- c(palette, colours[x])

      symbol <- cl%%25 ## There are only 25 diverse
      if(symbol == 4){
        symbol <- 0
      }
      symbols <- c(symbols, symbol)
    }
  }
  mapping <- list('palette' = palette, 'symbols'=symbols)
  return(mapping)
}
