df <- read.delim("~/evol-feat-snakemake/workflow/output/formatted_orthology_table.tsv")

CN_variance <- read.delim("~/evol-feat-snakemake/workflow/input/CN_variance.tab", header=FALSE)

print(CN_variance[tail(sort(CN_variance$V2)),])

topVariance <- head(CN_variance[order(CN_variance$V2, decreasing = T),])
topVariance

for(og in topVariance$V1){
  # og = '41813at6656'
  sdf <- df[df$orthogroup==og,]
  summary(sdf$species)
  hist(summary(sdf$species), breaks=100, main=og, xlab='per species copy-number')
}

sorted_CN_variance <- CN_variance[order(CN_variance$V2, decreasing = T),]
row.names(sorted_CN_variance) <- NULL
topVarianceOGs <- as.data.frame(sorted_CN_variance$V1)

write.table(topVarianceOGs, file = "~/Downloads/topVarianceOGs.txt", col.names = F, row.names = F, quote = F)
