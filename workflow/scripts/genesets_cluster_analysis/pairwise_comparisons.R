input_features <- na.omit(read.delim(snakemake@input[[1]]))

library(GGally, quietly=TRUE)

# ## Splitting dataframe into partitions of 4 variables as ggpairs plot is better
# ## visualized with few variables only at a time.

# set_membership_col <- input_features[, 2]
# df <- input_features[, 3:ncol(input_features)]

# n_variables_per_plot <- 4
# modulo <- ncol(df) %% n_variables_per_plot
# n_full_partitions <- (ncol(df) - modulo) / n_variables_per_plot
# n_variables_in_modulo_partition <- modulo
#
# partitions <- list()
# if(n_full_partitions > 0){
#   for(i in 1:n_full_partitions){
#     start = n_variables_per_plot * (i - 1) + 1
#     end = n_variables_per_plot * i
#     new_partition <- df[, start:end]
#     partitions[[i]] <- cbind.data.frame(set_membership_col, new_partition)
#   }
#   if(modulo > 0){
#     start = end + 1
#     end = end + modulo
#     modulo_partition <- df[, start:end]
#     partitions[[n_full_partitions + 1]] <- cbind.data.frame(set_membership_col, modulo_partition)
#   }
# } else {
#   start = 1
#   end = modulo
#   modulo_partition <- df[, start:end]
#   partitions[[1]] <- cbind.data.frame(set_membership_col, modulo_partition)
# }

sets <- input_features[, 2]
metrics <- input_features[, 3:ncol(input_features)]
scaled_metrics <- scale(metrics, center = T, scale = T)
df <- data.frame(cbind.data.frame(sets, metrics))
sdf <- data.frame(cbind.data.frame(sets, scaled_metrics))

print(sdf)

pdf(file=snakemake@output[[1]], height=30, width=30)
  ggpairs(df, aes(color=sets, alpha=0.4), cardinality_threshold = NULL)
dev.off()

pdf(file=snakemake@output[[2]], height=30, width=30)
  ggpairs(sdf, aes(color=sets, alpha=0.4), cardinality_threshold = NULL)
dev.off()
# for(i in 1:length(partitions)){
#   ggpairs(partitions[[i]], aes(color=set, alpha=0.4), cardinality_threshold = NULL)
# }
