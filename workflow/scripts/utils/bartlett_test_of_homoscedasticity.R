df <- merged_orthogroup_features[,2:ncol(merged_orthogroup_features)]

bartlett.test(df)

sdf <- scale(df)

for(i in 1:ncol(df)){
  print(colnames(df)[i])
  print(var(sdf[,i]))
}

bartlett.test(as.data.frame(sdf))
