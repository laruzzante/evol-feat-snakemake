df <- read.delim(snakemake@input[[1]])

source('utils.R')

pdf(file = snakemake@output[[1]])
plot_dendrogram_heatmap(df, aggregator_function='mean')
dev.off()
