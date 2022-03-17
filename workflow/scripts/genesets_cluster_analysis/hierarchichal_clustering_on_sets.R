source('scripts/cluster_analysis_functions.R')

input_features <- na.omit(read.delim(snakemake@input[[1]]))

cat('\nCounts table of orthogroups per set:\n')
print(table(input_features$set))
cat('\n')

scaled_features <- scale_features(input_features)
mean_matrix <- try(format_features_matrix(scaled_features, aggregator_function='mean'))
median_matrix <- try(format_features_matrix(scaled_features, aggregator_function='median'))

n_bootstraps <- 10000


pdf(file=snakemake@output[[1]])
  if(exists("mean_matrix")) {
    try(euclidean_ward.d2_mean <- plot_pvclust_basic_heatmap_on_sets(mean_matrix, nboot=n_bootstraps, aggregator_function='mean', method.dist='euclidean', method.hclust='ward.D2'))
    try(pearson_average_mean <- plot_pvclust_basic_heatmap_on_sets(mean_matrix, nboot=n_bootstraps, aggregator_function='mean', method.dist=pearson, method.hclust='average'))
  }
  if(exists("median_matrix")) {
    try(euclidean_ward.d2_median <- plot_pvclust_basic_heatmap_on_sets(median_matrix, nboot=n_bootstraps, aggregator_function='median', method.dist='euclidean', method.hclust='ward.D2'))
    try(pearson_average_median <- plot_pvclust_basic_heatmap_on_sets(median_matrix, nboot=n_bootstraps, aggregator_function='median', method.dist=pearson, method.hclust='average'))
  }
dev.off()
