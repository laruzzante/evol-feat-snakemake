library(Rtsne)
library(scales)
library(RColorBrewer)
library(uwot)

## Curating the dataset for analysis with both t-SNE and UMAP

set.seed(12345)
df <- na.omit(read.delim(snakemake@input[[1]]))
sdf <- scale(df[,-c(1,2)], center=T, scale=T)
train <- cbind(df[,c(1,2)], sdf)

## Definine colour palette for plots
N <- length(unique(train$set))
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
colours <- sample(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:N])
names(colours) <- unique(train$set)

## Executing the algorithm on curated data
perpl <- sqrt(nrow(sdf))
tsne <- Rtsne(sdf, dims = 2, perplexity=perpl, verbose=TRUE, max_iter=500, check_duplicates=FALSE)

umap <- try(umap(sdf, n_components = 2, n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8, metric = 'correlation'))
## If using correlation metric fails, then try with Euclidean distance
if(exists("umap") == FALSE){
  umap <- umap(sdf, n_components = 2, n_neighbors = 15, min_dist = 0.001, verbose = TRUE, n_threads = 8, metric = 'euclidean')
  print('WARNING: UMAP computation switched to Euclidean distance metric instead of Pearson Correlation due to computational errors.')
}

## Assigning colours and symbols to sets
unique_sets <- unique(train$set)
palette <- c()
symbols <- c()
for(set in train$set){
  x <- which(unique_sets == set)
  while(x > length(colours)){
    x <- x - length(colours)
  }
  palette <- c(palette, colours[x])
  symbol <- x%%25 ## Only 25 unique symbols in base-R pch
  symbols <- c(symbols, symbol)
}

## Plotting
pdf(snakemake@output[[1]])

plot(tsne$Y, t='n', main="tSNE on scaled metrics", xlab='tSNE.X', ylab='tSNE.Y')
points(tsne$Y, col=palette, cex=0.7, pch=symbols)
legend(x='topleft', legend=names(colours), pch=unique(symbols), col=colours, cex=0.7)

plot(umap, t='n', main="UMAP on scaled metrics")
points(umap, col=palette, cex=0.7, pch=symbols)
legend(x='topleft', legend=names(colours), pch=unique(symbols), col=colours, cex=0.7)

dev.off()
