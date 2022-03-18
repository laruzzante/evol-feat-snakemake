library(ggpubr)

set.seed(12345)

THREADS <- snakemake@config[['MAX_THREADS']]

merged_orthogroup_features <- read.delim(snakemake@input[['features']])
df <- na.omit(merged_orthogroup_features[,-1])
sdf <- scale(df, center=T, scale=T)

pdf(snakemake@output[['plot']], paper="a4")
plot.new()
txt.size <- 0.8
text.xpos <- 0.5
# Bartlett homoscedasticity test (homogeneity of variances)
bartlett <- bartlett.test(as.data.frame(sdf))
bartlett.txt <- paste('Scaled Metrics\n\n\n',bartlett$method,':\np-value: ',bartlett$p.value, sep='')
text(x=text.xpos, y=0.5, bartlett.txt, cex=txt.size)
for(i in 1:ncol(sdf)){
  hist(sdf[,i],breaks = 100, main='Scaled metric histogram', xlab=colnames(sdf)[i])
  # quantile-quantile plot
  print(ggqqplot(sdf[,i]))
}
dev.off()
