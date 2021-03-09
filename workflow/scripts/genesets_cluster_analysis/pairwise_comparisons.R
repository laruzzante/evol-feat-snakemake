input_features <- na.omit(read.delim(snakemake@input[[1]]))

library(GGally, quietly=TRUE)

pdf(file=snakemake@output[[1]], height=30, width=30)
  ggpairs(input_features[2:ncol(input_features)], aes(color=set, alpha=0.4))
dev.off()
