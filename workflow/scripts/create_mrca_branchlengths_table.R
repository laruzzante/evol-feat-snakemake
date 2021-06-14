## Building a lookup table to find out branch lengths of internal tree nodes
## in order to compute the AGE metric.

library(ape)

get_nodelengths_table <- function(tree){

  AGE_matrix <- matrix(rep(NA), Ntip(tree),Ntip(tree))
  rownames(AGE_matrix) <- tree$tip.label
  colnames(AGE_matrix) <- tree$tip.label

  AGE_dataframe <- c()

  for(spec1 in sort(tree$tip.label)){
    for(spec2 in sort(tree$tip.label)){
      if(spec1 != spec2){
        mrca_node <- as.numeric(getMRCA(phy=tree, tip = c(spec1, spec2)))
        mrca_age <- as.numeric(branching.times(tree)[mrca_node - Ntip(tree)])
        AGE_row <- c(spec1, spec2, mrca_age)
        AGE_dataframe <- rbind(AGE_dataframe, AGE_row)
      }
    }
  }
  colnames(AGE_dataframe) <- c("species1", "species2", "AGE")
  rownames(AGE_dataframe) <- 1:nrow(AGE_dataframe)
  return(as.data.frame(AGE_dataframe))
}

tree <- read.tree(snakemake@input[[1]])
table <- get_nodelengths_table(tree)

write.table(table, file=snakemake@output[[1]],
            sep='\t', quote=FALSE, row.name=FALSE)
