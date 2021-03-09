input_features <- na.omit(read.delim(snakemake@input[[1]]))

library(Rtsne)

train <- input_features

## Curating the database for analysis with both t-SNE and PCA
train$set<-as.factor(train$set)
## for plotting
colors = rainbow(length(unique(train$set)))
names(colors) = unique(train$set)

## Executing the algorithm on curated data
set.seed(1234)
tsne <- Rtsne(train[,-c(1,2)], dims = 2, perplexity=30, verbose=TRUE, max_iter=500, check_duplicates=FALSE)

## Plotting
pdf(snakemake@output[[1]])
  plot(tsne$Y, t='n', main="tsne")
  points(tsne$Y, col=colors[train$set], cex=0.8)
  legend(x='topright', legend=names(colors), pch=1, col=colors, cex=0.6)
dev.off()
