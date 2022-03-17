library(colorRamps)
library(dendextend)
library(gplots)
library(pvclust)

plot_pvclust_annotated_heatmap <- function(matrix, nboot=10000, aggregator_function='mean', method.dist='euclidean', method.hclust='ward.d2', iseed=1234){
  sets_dendrogram <- pvclust(t(matrix), method.dist=method.dist, method.hclust=method.hclust, nboot=nboot, parallel=T, iseed=iseed)
  features_dendrogram <- pvclust(matrix, method.dist=method.dist, method.hclust=method.hclust, nboot=nboot, parallel=T, iseed=iseed)

  sets <- as.dendrogram(sets_dendrogram$hclust)
  features <- as.dendrogram(features_dendrogram$hclust)

  # AU p-values from pvclust
  sets_AU <- sets_dendrogram$edges$au
  features_AU <- features_dendrogram$edges$au


  # Ordering the dendrogram nodes from pvclust following as.dendrogram node order, WHICH IS DIFFERENT!@#$%&*!!
  # And as.dendrogram also includes leaves indices within nodes. So be careful that leaves and nodes are all mixed togheter.
  # I recognize leaves by their height of 0.
  # I'm using height of nodes as the mapping variable from one pvclust dendrogram to as.dendrogram becuase it's the only
  # common factor I could use. So, pvclust dendrogram nodes are named by their height score and then ordered following the
  # height vector from as.dendrogram which corresponds to the actual plotting order of the nodes. The limitation is that
  # we must be careful in the rare eventuality that two nodes have the exact same height.
  names(sets_AU) <- sets_dendrogram$hclust$height # same for families
  dendrogram_sets_heights <- get_nodes_attr(sets, "height", include_leaves = T)
  sets_AU_sorted <- c()
  for(node_height in dendrogram_sets_heights){
    if(node_height != 0.0) sets_AU_sorted <- c(sets_AU_sorted, sets_AU[as.character(node_height)])
    else sets_AU_sorted <- c(sets_AU_sorted, NA)
  }

  names(features_AU) <- features_dendrogram$hclust$height # metrics AU p-values
  dendrogram_features_heights <- get_nodes_attr(features, "height", include_leaves = T)
  features_AU_sorted <- c()
  for(node_height in dendrogram_features_heights){
    if(node_height != 0.0) features_AU_sorted <- c(features_AU_sorted, features_AU[as.character(node_height)])
    else features_AU_sorted <- c(features_AU_sorted, NA)
  }

  ## Heatmap

  # Setting AU p-values colour-scale in dendrograms
  AU_colorpalette <- c(colorRampPalette(brewer.pal(9, "Greens")[3:9])(9))

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

  sets_colors_AU <- map2color(sets_AU_sorted)
  features_colors_AU <- map2color(features_AU_sorted)

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


  Rowv <- sets %>%
    #set("branches_k_color", k = 10) %>% set("branches_lwd", 3) %>%  # setting the number of separate colour-coded clades
    set("nodes_pch", 20) %>%  # node point type
    set("nodes_cex", 1.5) %>%  # node point size
    set("nodes_col", sets_colors_AU, warn=T) %>% # node point color
    set("clear_leaves") #%>%

  Colv <- features %>%
    #set("branches_k_color", k = 5) %>% set("branches_lwd", 3) %>%  # setting the number of separate colour-coded clades
    set("nodes_pch", 20) %>%  # node point type
    set("nodes_cex", 1.5) %>%  # node point size
    set("nodes_col", features_colors_AU, warn=T) %>% # node point color
    set("clear_leaves") #%>%

  hm <- heatmap.2(matrix, Rowv = Rowv, Colv = Colv, trace = "none", labRow = FALSE,
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

  return(hm)


input_features <- na.omit(read.delim(snakemake@input[[1]]))

source('scripts/genesets_cluster_analysis/functions.R')


cat('\nCounts table of orthogroups per set:\n')
print(table(input_features$set))
cat('\n')

scaled_features <- scale_features(input_features)
mean_matrix <- try(format_features_matrix(scaled_features, aggregator_function='mean'))
median_matrix <- try(format_features_matrix(scaled_features, aggregator_function='median'))

n_bootstraps <- 10000

pdf(file=snakemake@output[[1]])
  if(exists("mean_matrix")) {
    try(euclidean_ward.d2_mean <- plot_pvclust_annotated_heatmap(mean_matrix, nboot=n_bootstraps, aggregator_function='mean', method.dist='euclidean', method.hclust='ward.D2'))
    try(pearson_average_mean <- plot_pvclust_annotated_heatmap(mean_matrix, nboot=n_bootstraps, aggregator_function='mean', method.dist=pearson, method.hclust='average'))
  }
  if(exists("median_matrix")) {
    try(euclidean_ward.d2_median <- plot_pvclust_basic_heatmap(median_matrix, nboot=n_bootstraps, aggregator_function='median', method.dist='euclidean', method.hclust='ward.D2'))
    try(pearson_average_median <- plot_pvclust_basic_heatmap(median_matrix, nboot=n_bootstraps, aggregator_function='median', method.dist=pearson, method.hclust='average'))
  }
dev.off()
