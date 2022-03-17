## Check for normality
library(ggpubr)
library(nortest)
library(M3C)
library(Rtsne)
library(uwot)
setwd("~/evol-feat-snakemake/workflow/output")

merged_orthogroup_features <- read.delim("~/evol-feat-snakemake/workflow/output/merged_orthogroup_features.tsv")

df <- na.omit(merged_orthogroup_features[,2:ncol(merged_orthogroup_features)])

sdf <- scale(df)

set.seed(1234)
### Homoscedasticity and normality tests on:


## Scaled Metrics

pdf("~/evol-feat-snakemake/workflow/scaled_metrics_homoscedasticity_and_normality_tests.pdf", paper="a4")
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
  plot.new()
  # Shapiro-Wilk normality test
  shapirowilk <- shapiro.test(sample(sdf[,i],5000)) # Shapiro test only takes up to 5k data points, hence I am sampling a random 5k from the metric
  shapirowilk.txt <- paste(shapirowilk$method,':\np-value: ',shapirowilk$p.value, sep='')
  text(x=text.xpos, y=1, shapirowilk.txt, cex=txt.size)
  # Anderson-Darling test for normality
  andersondarling <- ad.test(sdf[,i])
  andersondarling.txt <- paste(andersondarling$method,':\np-value: ',andersondarling$p.value, sep='')
  text(x=text.xpos, y=0.7, andersondarling.txt, cex=txt.size)
  # KS test to check data's significant difference to any continuous density function 
  kolmogorowsmirnow <- ks.test(sdf[,i], y='pnorm')
  kolmogorowsmirnow.txt <- paste(kolmogorowsmirnow$method,':\np-value: ',kolmogorowsmirnow$p.value, sep='')
  text(x=text.xpos, y=0.4, kolmogorowsmirnow.txt, cex=txt.size)
}
dev.off()


## Principal Components Scores of Scaled Metrics

pc <- princomp(sdf)
pcdf <- pc$scores

pdf("~/evol-feat-snakemake/workflow/principal_components_scaled_metrics_homoscedasticity_and_normality_tests.pdf", paper="a4")
plot.new()
txt.size <- 0.8
text.xpos <- 0.5
# Bartlett homoscedasticity test (homogeneity of variances)
bartlett <- bartlett.test(as.data.frame(pcdf))
bartlett.txt <- paste('Principal Components Scores of Scaled Metrics\n\n\n',bartlett$method,':\np-value: ',bartlett$p.value, sep='')
text(x=text.xpos, y=0.5, bartlett.txt, cex=txt.size)
for(i in 1:ncol(pcdf)){
  hist(pcdf[,i],breaks = 100, main='PC of Scaled metric histogram', xlab=colnames(pcdf)[i])
  # quantile-quantile plot
  print(ggqqplot(pcdf[,i]))
  plot.new()
  # Shapiro-Wilk normality test
  shapirowilk <- shapiro.test(sample(pcdf[,i],5000)) # Shapiro test only takes up to 5k data points, hence I am sampling a random 5k from the metric
  shapirowilk.txt <- paste(shapirowilk$method,':\np-value: ',shapirowilk$p.value, sep='')
  text(x=text.xpos, y=1, shapirowilk.txt, cex=txt.size)
  # Anderson-Darling test for normality
  andersondarling <- ad.test(pcdf[,i])
  andersondarling.txt <- paste(andersondarling$method,':\np-value: ',andersondarling$p.value, sep='')
  text(x=text.xpos, y=0.7, andersondarling.txt, cex=txt.size)
  # KS test to check data's significant difference to any continuous density function 
  kolmogorowsmirnow <- ks.test(pcdf[,i], y='pnorm')
  kolmogorowsmirnow.txt <- paste(kolmogorowsmirnow$method,':\np-value: ',kolmogorowsmirnow$p.value, sep='')
  text(x=text.xpos, y=0.4, kolmogorowsmirnow.txt, cex=txt.size)
}
dev.off()


## Log of Principal Components Scores of Scaled Metrics

logdf <- log(pcdf+2*abs(min(pcdf)))

pdf("~/evol-feat-snakemake/workflow/log_principal_components_scaled_metrics_homoscedasticity_and_normality_tests.pdf", paper="a4")
plot.new()
txt.size <- 0.8
text.xpos <- 0.5
# Bartlett homoscedasticity test (homogeneity of variances)
bartlett <- bartlett.test(as.data.frame(logdf))
bartlett.txt <- paste('Log Principal Components Scores of Scaled Metrics\n\n\n',bartlett$method,':\np-value: ',bartlett$p.value, sep='')
text(x=text.xpos, y=0.5, bartlett.txt, cex=txt.size)
for(i in 1:ncol(logdf)){
  hist(logdf[,i],breaks = 100, main='PC of Scaled metric histogram', xlab=colnames(logdf)[i])
  # quantile-quantile plot
  print(ggqqplot(logdf[,i]))
  plot.new()
  # Shapiro-Wilk normality test
  shapirowilk <- shapiro.test(sample(logdf[,i],5000)) # Shapiro test only takes up to 5k data points, hence I am sampling a random 5k from the metric
  shapirowilk.txt <- paste(shapirowilk$method,':\np-value: ',shapirowilk$p.value, sep='')
  text(x=text.xpos, y=1, shapirowilk.txt, cex=txt.size)
  # Anderson-Darling test for normality
  andersondarling <- ad.test(logdf[,i])
  andersondarling.txt <- paste(andersondarling$method,':\np-value: ',andersondarling$p.value, sep='')
  text(x=text.xpos, y=0.7, andersondarling.txt, cex=txt.size)
  # KS test to check data's significant difference to any continuous density function 
  kolmogorowsmirnow <- ks.test(logdf[,i], y='pnorm')
  kolmogorowsmirnow.txt <- paste(kolmogorowsmirnow$method,':\np-value: ',kolmogorowsmirnow$p.value, sep='')
  text(x=text.xpos, y=0.4, kolmogorowsmirnow.txt, cex=txt.size)
}
dev.off()


# ## tSNE Scores of Scaled Metrics
# 
# perplexity <- floor((nrow(sdf) - 1) / 3)
# tsnedf <- Rtsne(sdf, check_duplicates=FALSE, pca = F, dims = 3, num_threads = 0) # Defaults: perplexity=30 and pca=TRUE ... i.e. we have to specify that we do not want a PCA beforehand
# 
# pdf("~/evol-feat-snakemake/workflow/log_principal_components_scaled_metrics_homoscedasticity_and_normality_tests.pdf", paper="a4")
# plot.new()
# txt.size <- 0.8
# text.xpos <- 0.5
# 
# # Bartlett homoscedasticity test (homogeneity of variances)
# bartlett <- bartlett.test(as.data.frame(logdf))
# bartlett.txt <- paste('Log Principal Components Scores of Scaled Metrics\n\n\n',bartlett$method,':\np-value: ',bartlett$p.value, sep='')
# text(x=text.xpos, y=0.5, bartlett.txt, cex=txt.size)
# for(i in 1:ncol(tsnedf)){
#   hist(logdf[,i],breaks = 100, main='PC of Scaled metric histogram', xlab=colnames(logdf)[i])
#   # quantile-quantile plot
#   print(ggqqplot(logdf[,i]))
#   plot.new()
#   # Shapiro-Wilk normality test
#   shapirowilk <- shapiro.test(sample(logdf[,i],5000)) # Shapiro test only takes up to 5k data points, hence I am sampling a random 5k from the metric
#   shapirowilk.txt <- paste(shapirowilk$method,':\np-value: ',shapirowilk$p.value, sep='')
#   text(x=text.xpos, y=1, shapirowilk.txt, cex=txt.size)
#   # Anderson-Darling test for normality
#   andersondarling <- ad.test(logdf[,i])
#   andersondarling.txt <- paste(andersondarling$method,':\np-value: ',andersondarling$p.value, sep='')
#   text(x=text.xpos, y=0.7, andersondarling.txt, cex=txt.size)
#   # KS test to check data's significant difference to any continuous density function 
#   kolmogorowsmirnow <- ks.test(logdf[,i], y='pnorm')
#   kolmogorowsmirnow.txt <- paste(kolmogorowsmirnow$method,':\np-value: ',kolmogorowsmirnow$p.value, sep='')
#   text(x=text.xpos, y=0.4, kolmogorowsmirnow.txt, cex=txt.size)
# }
# dev.off()
