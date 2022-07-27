library(GGally)
library(cowplot)

merged_orthogroup_features <- read.delim("//wsl.localhost/Ubuntu-20.04/home/lruzzant/evol-feat-snakemake/workflow/output/merged_orthogroup_features.tsv")
df <- na.omit(merged_orthogroup_features[,-1])
sdf <- as.data.frame(scale(df))

ggpairs(sdf, lower=list(continuous=wrap("points",alpha=0.1,size=0.2,col='blue')),
             diag=list(continuous=wrap("densityDiag",col='red',lwd=1)),
             upper=list(continuous=wrap("cor",size=4,col='black'))) +
    theme_bw(base_size = 8)

# ggsave("pairwise_comparisons.png")
