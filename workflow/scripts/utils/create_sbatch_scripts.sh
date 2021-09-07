#!/usr/bin/env bash

while read og;
do
	grep -vw "${og}" test_topVarianceOGs.txt > "filter_lists/${og}_topVarianceOGs.txt"
	grep -vwf "filter_lists/${og}_topVarianceOGs.txt" arthropoda_gene_counts.tsv > "gene_count_tables/${og}_arthropoda_gene_counts.tsv"
 	sed "s/orthogroup/${og}/g" sbatch_def.sh > sbatch_scripts/sbatch_"$og".sh 
done < test_topVarianceOGs.txt

