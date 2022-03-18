
source('scripts/cluster_analysis_functions.R')

library(dbscan)

tsne.coords <- read.delim(snakemake@input[['tsne']], header = F)
merged_orthogroup_features <- read.delim(snakemake@input[['features']])

minPoints <- 5

## OPTICS on tSNE
pdf(file=snakemake@output[['plot']])

dbscan::kNNdistplot(tsne.coords, k = minPoints)
abline(h = 0.3, lty = 2) # The knee seems to be around 0.3, hence we plot a line at h=0.3 just to see the actual intersection
cl.optics <- optics(tsne.coords, eps = 0.3, minPts = minPoints)
cl.optics.cut <- extractDBSCAN(cl.optics, eps_cl = 0.3)
mapping <- map_palette_to_clusters(cl.optics.cut)

plot(tsne.coords, col=mapping$palette, pch=mapping$symbols, cex=0.2, lwd=0.2,
     xlab='tSNE.X', ylab='tSNE.Y')
dev.off()

orthogroups <- na.omit(merged_orthogroup_features)[,1]
cluster_ids <- cl.optics.cut$cluster
memberships <- cbind(orthogroups, cluster_ids)

write.table(memberships, snakemake@output[['optics_clusters']], quote = F, row.names = F, col.names = F, sep='\t')
