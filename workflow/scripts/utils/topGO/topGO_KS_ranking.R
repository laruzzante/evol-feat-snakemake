set.seed(12345)

setwd("//wsl.localhost/Ubuntu-20.04/home/lruzzant/evol-feat-snakemake/")
output_dir <- "//wsl.localhost/Ubuntu-20.04/home/lruzzant/evol-feat-data-processing/topGO_ks/" 


merged_orthogroup_features <- read.delim("//wsl.localhost/Ubuntu-20.04/home/lruzzant/evol-feat-snakemake/workflow/output/merged_orthogroup_features.tsv")
df <- na.omit(merged_orthogroup_features[,-1])

library(topGO, quietly = T)
library(Rgraphviz, quietly = T)

go_universe <- 'workflow/output/go_enrichment_analysis/orthogroups_go_universe.tsv'
ontology <- 'BP'

print(ontology)

orthogroups2go <- readMappings(file = go_universe)
orthogroupsUniverse <- names(orthogroups2go)
str(head(orthogroups2go))

for(sel in c('topQuantile', 'bottomQuantile')){
  print(sel)
  
  for(metric in colnames(df)){
    print(metric)
    
    ogList <- merged_orthogroup_features[,metric]
    names(ogList) <- na.omit(merged_orthogroup_features[,"orthogroup"])
    ogList <- na.omit(ogList)
    
    # Gene selection function
    topcut <- 0.9
    
    mySelgn <- function(score) {
      if(sel == 'min') return(score <= (max(ogList)*(1-topcut)) )
      else if(sel == 'max') return(score >= (max(ogList)*topcut) )
      else if(sel =='bottomQuantile') return(score <= quantile(ogList, 1-topcut, names=FALSE))
      else if(sel == 'topQuantile') return(score >= quantile(ogList, topcut, names=FALSE))
    }
    
    GOdata <- new("topGOdata", ontology = ontology, allGenes = ogList,
                  annot = annFUN.gene2GO, gene2GO = orthogroups2go, geneSel = mySelgn, nodeSize = 10)
    
    
    
    if(sel == "topQuantile" | sel == "max"){
      scoreOrder <- "decreasing"
    } else if(sel =="bottomQuantile" | sel == "min"){
      scoreOrder <- "inreasing"
    }
    
    resKSweight01 <- runTest(GOdata, algorithm = "weight01", statistic = "ks", scoreOrder = scoreOrder)
    print(resKSweight01)
    
    weight01Res <- GenTable(GOdata, W01KS = resKSweight01, ranksOf = "W01KS", topNodes = 50)
    print(weight01Res)
    
    output_table <- paste0(output_dir, metric, '_', as.character(1-topcut), sel, '.tsv')
    
    write.table(weight01Res, file=output_table, quote = F, sep='\t', col.names = T, row.names = F, append=T)
  }
}
