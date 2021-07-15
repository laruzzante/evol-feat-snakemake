df <- read.delim("~/evol-feat-snakemake/workflow/output/merged_orthogroup_features.tsv")

library(plot3D)

i <- 2
j <- 8

x <- df[[i]]
y <- df[[j]]

##  Create cuts:
n_cuts <- 50
x_c <- cut(x, n_cuts)
y_c <- cut(y, n_cuts)

##  Calculate joint counts at cut levels:
z <- table(x_c, y_c)


xlab <- paste0(colnames(df)[i], ': [',as.character(round(min(df[[i]]),2)), ', ', as.character(round(max(df[[i]]),2)), ']')
ylab <- paste0(colnames(df)[j], ': [',as.character(round(min(df[[j]]),2)), ', ', as.character(round(max(df[[j]]),2)), ']')

##  Plot as a 2D heatmap:
image2D(z=z, border="black", xlab=xlab, ylab=ylab)

##  Plot as a 3D histogram:
hist3D(z=z, border='black', xlab=xlab, ylab=ylab)
