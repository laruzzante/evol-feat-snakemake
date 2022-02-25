library(pvclust)
library(ComplexHeatmap)
source('~/agambiae-immunity-families/livio/scripts/functions.R')

n_sample <- 200

merged_orthogroup_features <- read.delim("~/Downloads/merged_orthogroup_features.tsv")

df <- merged_orthogroup_features[,2:ncol(merged_orthogroup_features)]

df_s <- df[sample(1:nrow(df), n_sample), ]

mat <- as.matrix(df_s)

evol_metrics_matrix <- scale(mat)


bootstraps_metrics <- pvclust(evol_metrics_matrix, method.dist = pearson, method.hclust = "complete", nboot = 100, parallel = T, store = F, iseed = 1234)
bootstraps_families <- pvclust(t(evol_metrics_matrix), method.dist = pearson, method.hclust = "complete", nboot = 100, parallel = T, store = F, iseed = 1234)


heatmap_colors <- rev(unlist(strsplit('1E009F-2531C7-0056FF-3C8DF9-6BAAF0-BCDBF7-F8F987-FFCA00-F9A500-FF6500-FF3800-D20000', '-')))
heatmap_colors <- paste0("#", heatmap_colors)

# hm <- heatmap.2(evol_metrics_matrix, Rowv = bootstraps_families$hclust, Colv = bootstraps_metrics$hclust, trace = "none",
#                 density.info = "none", col = rev(heatmap_colors), scale = "column", key = FALSE)

library(circlize)
# col_fun = colorRamp2(seq(-2, 2), rev(c("red", "yellow", "white", "cyan", "darkblue")))
col_fun = colorRamp2(seq(-5.5, 5.5), rev(heatmap_colors))

plot(bootstraps_metrics)
plot(bootstraps_families)
Heatmap(evol_metrics_matrix, cluster_rows = bootstraps_families$hclust, cluster_columns = bootstraps_metrics$hclust, 
        col = col_fun, heatmap_legend_param = list(title = "Scaled\nMetrics"), show_row_names = F)

# hm <- heatmap.2(evol_metrics_matrix, Rowv = Rowv, Colv = Colv, trace = "none",
#                 density.info = "none", col = rev(heatmap_colors), scale = "column", key = FALSE,
#                 main = paste0(aggregator_function, '\n', 'method.dist: ', distance_method, '\nmethod.hclust: ', hclust_method))
# 
# hm <- Heatmap(evol_metrics_matrix, cluster_rows = bootstraps_families$hclust, cluster_columns = bootstraps_metrics$hclust, heatmap_legend_param = list(title = "Scaled\nMetrics"))
