input_features <- na.omit(read.delim(snakemake@input[[1]]))

source('scripts/cluster_analysis_functions.R')
library(factoextra, quietly=TRUE)

scaled_features <- scale_features(input_features)
pc_orthogroups <- prcomp(scaled_features[3:ncol(scaled_features)], center=TRUE, scale.=TRUE)


n_principal_components <- 10
n_bootstraps <- 10000
seed <- 1234
method.dist <- 'euclidean'
method.hclust <- 'ward.D2'

mean_matrix <- format_features_matrix(scaled_features, aggregator_function='mean')
pc_mean_sets <- prcomp(mean_matrix, center=TRUE, scale.=TRUE)$x[,1:n_principal_components]
mean_sets_dendrogram <- pvclust(t(pc_mean_sets), method.dist=method.dist, method.hclust=method.hclust, nboot=n_bootstraps, parallel=T, iseed=seed)

# median_matrix <- format_features_matrix(scaled_features, aggregator_function='median')
# pc_median_sets <- prcomp(median_matrix, center=TRUE, scale.=TRUE)$x[,1:n_principal_components]
# median_sets_dendrogram <- pvclust(t(pc_median_sets), method.dist=method.dist, method.hclust=method.hclust, nboot=n_bootstraps, parallel=T, iseed=seed)


pdf(file=snakemake@output[[1]])

  # Plotting PCA for all orthogroups, and barplot of the first two principal components
  fviz_pca_biplot(pc_orthogroups, label = 'var', habillage=scaled_features$set, col.var='#ff5e5e', addEllipses=TRUE)
  barplot(rev(sort(pc_orthogroups$rotation[,1])), main='PC1', las=2)
  barplot(rev(sort(pc_orthogroups$rotation[,2])), main='PC2', las=2)

  # Plotting pvclust dendrograms of first N principal components
  plot(mean_sets_dendrogram, main=paste('orthogroup sets averaged by mean\n',
                                        n_principal_components, 'principal components'))
  pvrect(mean_sets_dendrogram)
  # plot(median_sets_dendrogram, main=paste('orthogroup sets averaged by median\n',
  #                                       n_principal_components, 'principal components'))
  # pvrect(median_sets_dendrogram)

dev.off()
