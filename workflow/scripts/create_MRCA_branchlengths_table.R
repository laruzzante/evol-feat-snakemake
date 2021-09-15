## Building a lookup table to find out branch lengths of internal tree nodes
## in order to compute the branchlength metric.

library(ape)

get_nodelengths_table <- function(tree){

  branchlength_matrix <- matrix(rep(NA), Ntip(tree),Ntip(tree))
  rownames(branchlength_matrix) <- tree$tip.label
  colnames(branchlength_matrix) <- tree$tip.label

  branchlength_dataframe <- c()

  for(spec1 in sort(tree$tip.label)){
    for(spec2 in sort(tree$tip.label)){
      if(spec1 != spec2){
        mrca_node <- as.numeric(getMRCA(phy=tree, tip = c(spec1, spec2)))
        mrca_branchlength <- as.numeric(branching.times(tree)[mrca_node - Ntip(tree)])
        branchlength_row <- c(spec1, spec2, mrca_branchlength)
        branchlength_dataframe <- rbind(branchlength_dataframe, branchlength_row)
      }
    }
  }
  colnames(branchlength_dataframe) <- c("species1", "species2", "branchlength")
  rownames(branchlength_dataframe) <- 1:nrow(branchlength_dataframe)
  return(as.data.frame(branchlength_dataframe))
}

tree <- read.tree(snakemake@input[[1]])
table <- get_nodelengths_table(tree)

write.table(table, file=snakemake@output[[1]],
            sep='\t', quote=FALSE, row.name=FALSE)
