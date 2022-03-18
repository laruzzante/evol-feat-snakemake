library(nortest)

set.seed(12345)

THREADS <- snakemake@config[['MAX_THREADS']]

merged_orthogroup_features <- read.delim(snakemake@input[['features']])
df <- na.omit(merged_orthogroup_features[,-1])
sdf <- scale(df, center=T, scale=T)

pdf(snakemake@output[['plot']], paper="a4")
for(i in 1:ncol(sdf)){
  plot.new()
  txt.size <- 0.8
  text.xpos <- 0.5
  text(x=text.xpos, y=1, colnames(sdf)[i], cex=txt.size)
  # Shapiro-Wilk normality test
  shapirowilk <- shapiro.test(sample(sdf[,i],5000)) # Shapiro test only takes up to 5k data points, hence I am sampling a random 5k from the metric
  shapirowilk.txt <- paste(shapirowilk$method,':\np-value: ',shapirowilk$p.value, sep='')
  text(x=text.xpos, y=0.8, shapirowilk.txt, cex=txt.size)
  # Anderson-Darling test for normality
  andersondarling <- ad.test(sdf[,i])
  andersondarling.txt <- paste(andersondarling$method,':\np-value: ',andersondarling$p.value, sep='')
  text(x=text.xpos, y=0.6, andersondarling.txt, cex=txt.size)
  # KS test to check data's significant difference to any continuous density function
  kolmogorowsmirnow <- ks.test(sdf[,i], y='pnorm')
  kolmogorowsmirnow.txt <- paste(kolmogorowsmirnow$method,':\np-value: ',kolmogorowsmirnow$p.value, sep='')
  text(x=text.xpos, y=0.4, kolmogorowsmirnow.txt, cex=txt.size)
}
dev.off()
