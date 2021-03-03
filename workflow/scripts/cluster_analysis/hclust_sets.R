df <- read.delim(snakemake@input[[1]])

aggregator_function = 'mean'

scale_metrics <- function(df){
  scaled_df <- df
  for(i in 1:ncol(df)){
    if(sapply(df[i], is.numeric)){
      scaled_df[i] <- scale(df[i])
    } else cat(paste('Warning/check:', colnames(df[i]), "-> non-numeric, not scalable\n"))
  }
  return(scaled_df)
}

format_metrics_matrix <- function(scaled_df, aggregator_function=aggregator_function){
  df <- scaled_df
  df <- aggregate(df[,3:ncol(df)], by=list(df[['set']]), FUN = aggregator_function, na.rm=TRUE)
  rownames(df) <- df$Group.1
  evol_metrics_matrix <- as.matrix(df[2:ncol(df)])
  return(evol_metrics_matrix)
}

mat <- scale_metrics(df)
mat <- format_metrics_matrix(df, aggregator_function = aggregator_function)

# bootstraps_metrics <- pvclust(df, method.dist = method.dist, method.hclust = method.hclust, nboot = 100, parallel = T, iseed = 1234)

install.packages("gplots", repos="https://cran.rstudio.com")
library(gplots)
install.packages("RColorBrewer", repos="https://cran.rstudio.com")
library(RColorBrewer)

# hclust(df)
#
# summary(df)
#
# plot(hclust(dist(mat), member=table(df$set)))
# plot(hclust(dist(mat)))
hmcol <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(150))

pdf(file = snakemake@output[[1]])

hm <- heatmap.2(mat, scale="column",  col=hmcol,
                main=paste('Features', aggregator_function), na.color='grey', trace='none' , margins=c(11,13), srtCol=40)

dev.off()
