#!/usr/bin/env Rscript

## R script which returns the branchlength of the LCA between two species in a phylogenetic tree

library(ape)

setwd('~/evol-feat/')

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least three arguments must be supplied (input file).
       Please specify a phylogenetic tree readable by ape and
       at least two OrthoDB species codenames present in the tree.n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

tree <- read.tree(args[1])

max_age <- 0
for(spec1 in args[2:length(args)]){
  for(spec2 in args[3:length(args)]){
    mrca_node <- as.numeric(getMRCA(phy=tree, tip = c(spec1, spec2)))
    mrca_age <- as.numeric(branching.times(tree)[mrca_node - Ntip(tree)])
    if(mrca_age > max_age){
      max_age <- mrca_age
    }
  }
}

print(max_age)
