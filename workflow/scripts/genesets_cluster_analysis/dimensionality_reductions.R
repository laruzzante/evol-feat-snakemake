input_features <- na.omit(read.delim(snakemake@input[[1]]))

library(Rtsne)
library(scales)
library(RColorBrewer)
library(uwot)

## Curating the database for analysis with both t-SNE and PCA

n_principal_components <- 10
pc_mean_sets <- prcomp(input_features[,-c(1,2)], center=TRUE, scale.=TRUE)$x[,1:n_principal_components]
train <- cbind(input_features[,c(1,2)], pc_mean_sets)
train$set<-as.factor(train$set)

## Definine colour palette for plots
N <- length(unique(train$set))
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
colors <- sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:N])
names(colors) <- unique(train$set)

## Executing the algorithm on curated data
set.seed(1234)
tsne_on_pc <- Rtsne(train[,-c(1,2)], dims = 2, perplexity=30, verbose=TRUE, max_iter=500, check_duplicates=FALSE)
tsne <- Rtsne(scale(input_features[,-c(1,2)]), dims = 2, perplexity=30, verbose=TRUE, max_iter=500, check_duplicates=FALSE)

umap_on_pc <- umap(train[,-c(1,2)], n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8)
umap <- umap(scale(input_features[,-c(1,2)]), n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8, metric = 'correlation')


## Plotting
pdf(snakemake@output[[1]])
plot(tsne_on_pc$Y, t='n', main="tSNE on first 10 PCs")
points(tsne_on_pc$Y, col=colors[tsne_on_pc$set], bg=alpha(colors[train$set],0.6), cex=1, pch=21)
legend(x='topleft', legend=names(colors), pch=21, col=colors, cex=0.4)

plot(umap_on_pc, t='n', main="UMAP on first 10 PCs")
points(umap_on_pc, col=colors[train$set], bg=alpha(colors[train$set],0.6), cex=1, pch=21)
legend(x='topleft', legend=names(colors), pch=21, col=colors, cex=0.4)

plot(tsne$Y, t='n', main="tSNE on scaled metrics")
points(tsne$Y, col=colors[train$set], bg=alpha(colors[train$set],0.6), cex=1, pch=21)
legend(x='topleft', legend=names(colors), pch=21, col=colors, cex=0.4)

plot(umap, t='n', main="UMAP on scaled metrics")
points(umap, col=colors[train$set], bg=alpha(colors[train$set],0.6), cex=1, pch=21)
legend(x='topleft', legend=names(colors), pch=21, col=colors, cex=0.4)
dev.off()
