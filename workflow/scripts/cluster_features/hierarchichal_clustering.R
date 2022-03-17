
source('scripts/cluster_analysis_functions.R')

merged_orthogroup_features <- read.delim(snakemake@input[[1]])
df <- na.omit(merged_orthogroup_features[,2:ncol(merged_orthogroup_features)])

## Scaling and centering
sdf <- scale(df, center = TRUE, scale = TRUE)

matrix <- as.matrix(sdf)
matrix_scaled_by_row <- t(scale(t(matrix)))
n_bootstraps <- 100

pdf(file=snakemake@output[[1]])
pearson_average <- plot_pvclust_basic_heatmap(matrix, nboot=n_bootstraps, scaling='columns', method.dist=pearson, method.hclust='average')
pearson_average <- plot_pvclust_basic_heatmap(matrix_scaled_by_row, scaling='columns and rows', nboot=n_bootstraps, method.dist=pearson, method.hclust='average')
dev.off()
