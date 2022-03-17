library(factoextra)

merged_orthogroup_features <- read.delim(snakemake@input[[1]])
df <- na.omit(merged_orthogroup_features[,2:ncol(merged_orthogroup_features)])

## Scaling and centering
sdf <- scale(df, center = TRUE, scale = TRUE)

## Principal Components to use as weighted scores for DBSCAN
pc <- princomp(sdf)

print(summary(pc))
print(pc$loadings)

## Printing PCA plots

pdf(file=snakemake@output[['plot']])

# Control variable colors using their contributions
# Control variable colors using their contributions
fviz_pca_var(pc, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
) + theme_minimal()

fviz_pca_biplot(pc, pointsize=0.2, label='var', repel=TRUE, col.var = 'red',
                col.ind='blue', alpha.ind=0.2) + theme_minimal()

# Contributions of variables
for(i in 1:length(pc$sdev)){
  print(fviz_contrib(pc, choice = "var", axes = i, top = Inf))
}

dev.off()
