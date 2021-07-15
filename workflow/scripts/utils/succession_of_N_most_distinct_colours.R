library(RColorBrewer)

n <- 50
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

colours <- col
pie(rep(1,n), col=colours)
