input_features <- na.omit(read.delim(snakemake@input[[1]]))

library(Rtsne)
library(scales)
library(RColorBrewer)

## Curating the database for analysis with both t-SNE and PCA

n_principal_components <- 5
pc_mean_sets <- prcomp(input_features[,-c(1,2)], center=TRUE, scale.=TRUE)$x[,1:n_principal_components]
train <- cbind(input_features[,c(1,2)],pc_mean_sets)
train$set<-as.factor(train$set)
## for plotting
N <- length(unique(train$set))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colors = sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:N])
names(colors) = unique(train$set)

## Executing the algorithm on curated data
set.seed(1234)
tsne <- Rtsne(train[,-c(1,2)], dims = 2, perplexity=30, verbose=TRUE, max_iter=500, check_duplicates=FALSE)

## Plotting
pdf(snakemake@output[[1]])
plot(tsne$Y, t='n', main="tsne")
points(tsne$Y, col=colors[train$set], bg=alpha(colors[train$set],0.6), cex=1, pch=21)
legend(x='topleft', legend=names(colors), pch=19, col=colors, cex=0.4)
dev.off()
