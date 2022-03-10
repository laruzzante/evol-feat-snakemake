# Succession of N most distinct colours
n <- 300
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
colours = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


cl.optics <- optics(tsnedf$Y, eps = 0.3)
cl.optics.cut <- extractDBSCAN(cl.optics, eps_cl = 0.3)
palette <- c()
for(cl in cl.optics.cut$cluster){
  x <- cl + 1L # So to never have the cluster 0, but 1, because we need to start the indexing at 1 to extract the first colour in colours
  while(x > length(colours)){
    x <- x - length(colours)
  }
  palette <- c(palette, colours[x])
}
plot(tsnedf$Y, pch=1, col=palette, cex=0.2)
