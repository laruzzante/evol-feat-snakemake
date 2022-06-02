source('C:/Users/livio/evol-feat-snakemake/workflow/scripts/cluster_analysis_functions.R')


setwd("C:/Users/livio/Desktop/ogsets/supercl")
df <- read.delim("C:/Users/livio/Desktop/ogsets/supercl/orthogroup_sets_features.tsv")

library(dplyr)

sets <- unique(sort(df$set))
features <- colnames(df)[-c(1,2)]

dat <- as.data.frame(matrix(ncol=length(features)+2, nrow=length(sets)))

colnames(dat) <- c('set', 'size', features)
dat$set <- sets

ssize <- df %>% group_by(set) %>% summarize(size = length(orthogroup))
dat$size <- ssize$size

for(feature in features){
  smeans <- c()
  for(set in sets){
    sset <- df[df$set==set,]
    smean <- mean(sset[,feature])
    smeans <- c(smeans, smean)
  }
  dat[,feature] <- smeans
}

write.table(dat, file='supercl_means.tsv', sep='\t', quote = F, row.names = F)
write.table(scale(dat), file='scaled_supercl_means.tsv', sep='\t', quote = F, row.names = F)


matrix <- scale(as.matrix(dat[,3:length(colnames(dat))]))
rownames(matrix) <- sets



input_features <- na.omit(read.delim(snakemake@input[[1]]))

cat('\nCounts table of orthogroups per set:\n')
print(table(input_features$set))
cat('\n')

scaled_features <- scale_features(input_features)
mean_matrix <- try(format_features_matrix(scaled_features, aggregator_function='mean'))
median_matrix <- try(format_features_matrix(scaled_features, aggregator_function='median'))

n_bootstraps <- 100


pdf(file=snakemake@output[[1]])
if(exists("mean_matrix")) {
  try(euclidean_ward.d2_mean <- plot_pvclust_basic_heatmap_on_sets(matrix, nboot=n_bootstraps, aggregator_function='mean', method.dist='euclidean', method.hclust='ward.D2'))
  try(pearson_average_mean <- plot_pvclust_basic_heatmap_on_sets(matrix, nboot=n_bootstraps, aggregator_function='mean', method.dist=pearson, method.hclust='average'))
}
if(exists("median_matrix")) {
  try(euclidean_ward.d2_median <- plot_pvclust_basic_heatmap_on_sets(median_matrix, nboot=n_bootstraps, aggregator_function='median', method.dist='euclidean', method.hclust='ward.D2'))
  try(pearson_average_median <- plot_pvclust_basic_heatmap_on_sets(median_matrix, nboot=n_bootstraps, aggregator_function='median', method.dist=pearson, method.hclust='average'))
}
dev.off()



scale_features <- function(df){
  scaled_df <- df
  for(i in 1:ncol(df)){
    if(sapply(df[i], is.numeric)){
      scaled_df[i] <- scale(df[i], center=TRUE, scale=TRUE)
    }
  }
  return(scaled_df)
}
