
source('scripts/cluster_analysis_functions.R')

library(dbscan)

umap.coords <- read.delim(snakemake@input[['umap']], header = F)
merged_orthogroup_features <- read.delim(snakemake@input[['features']])

minPoints <- 10

## OPTICS on tSNE
pdf(file=snakemake@output[['plot']])

dbscan::kNNdistplot(umap.coords, k = minPoints)
abline(h = 0.15, lty = 2) # The knee seems to be around 0.3, hence we plot a line at h=0.3 just to see the actual intersection
cl.hdbscan <- hdbscan(umap.coords, minPts = minPoints, gen_hdbscan_tree = T)

mapping <- map_palette_to_clusters(hdbscan)

plot(umap.coords, col=mapping$palette, pch=mapping$symbols, cex=0.2, lwd=0.2)
dev.off()

orthogroups <- na.omit(merged_orthogroup_features)[,1]
cluster_ids <- cl.hdbscan$cluster
memberships <- cbind(orthogroups, cluster_ids)

write.table(memberships, snakemake@output[['hdbscan_clusters']], quote = F, row.names = F, col.names = F, sep='\t')
