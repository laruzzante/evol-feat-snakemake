## Script that computes the number of children (tips/leaves) from any Most Recent Common Ancestor (internal nodes) in the phylogeny.

library(ape)

get_mrca_Ntip <- function(tree){

  subtrees <- subtrees(tree)

  Ntip_matrix <- matrix(rep(NA), Ntip(tree),Ntip(tree))
  rownames(Ntip_matrix) <- tree$tip.label
  colnames(Ntip_matrix) <- tree$tip.label

  Ntip_dataframe <- c()

  for(spec1 in sort(tree$tip.label)){
    for(spec2 in sort(tree$tip.label)){
      if(spec1 != spec2){
        mrca_index <- as.numeric(getMRCA(phy=tree, tip = c(spec1, spec2)))

        # Tips are labeled from [1] to [Number_of_leaves]. Nodes are labeled from
        # [N_number_of_leaves + 1] to [Number_of_nodes]. However nodes in the subtrees are
        # indexed starting from 1, hence we need to substract the total number of Tips
        # from the MRCA_index in order to correctly access the equivalent subtree from
        # the subtrees function.
        subtree_index = mrca_index - length(tree$tip.label)
        mrca_subtree <- subtrees[[subtree_index]]
        Ntip <- mrca_subtree$Ntip
        Ntip_row <- c(spec1, spec2, Ntip)
        Ntip_dataframe <- rbind(Ntip_dataframe, Ntip_row)
      }
    }
  }
  colnames(Ntip_dataframe) <- c("species1", "species2", "Ntip")
  rownames(Ntip_dataframe) <- 1:nrow(Ntip_dataframe)
  return(as.data.frame(Ntip_dataframe))
}

tree <- read.tree(snakemake@input[[1]])
table <- get_mrca_Ntip(tree)

write.table(table, file=snakemake@output[[1]],
            sep='\t', quote=FALSE, row.name=FALSE)
