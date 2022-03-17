
source('scripts/cluster_analysis_functions.R')

library(dbscan)

merged_orthogroup_features <- read.delim(snakemake@input[[1]])
df <- na.omit(merged_orthogroup_features[,2:ncol(merged_orthogroup_features)])

minPoints <- 10

## Scaling and centering
sdf <- scale(df, center = TRUE, scale = TRUE)

## Principal Components to use as weighted scores for DBSCAN
pc <- princomp(sdf)

## DBSCAN on weighted Principal Components
pdf(file=snakemake@output[['plot']])
dbscan::kNNdistplot(pc$scores, k = minPoints) ## By plotting this one we can check where the knee is. Careful, k must be equal to minPoints used above in dbscan
abline(h = 1.5, lty = 2) # The knee seems to be around 2, hence we plot a line at h=2 just to see the actual intersection
cl.dbscan <- dbscan(pc$scores, eps=1.5, weights = pc$sdev.^2 / sum(pc$sdev), minPts = minPoints)
# I have added weights on the clustering, i.e. each PC is being weighted by the proportion of explained variance
mapping <- map_palette_to_clusters(cl.dbscan)

plot(pc$scores[,c(1,2)], col=mapping$palette, pch=mapping$symbols, cex=0.2, lwd=0.2)
dev.off()

orthogroups <- na.omit(merged_orthogroup_features)[,1]
cluster_ids <- cl.dbscan$cluster
memberships <- cbind(orthogroups, cluster_ids)

write.table(memberships, snakemake@output[['weighted_dbscan_clusters']], quote = F, row.names = F, col.names = F, sep='\t')
