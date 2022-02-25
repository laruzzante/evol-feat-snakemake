21.11.2019 - Livio Ruzzante

1st step project: Arthur Royston, Waterhouse group

1) Initial trees: fungi.tre; insects.tre; nematodes.tre (i.e. *.tre)

2) Rooted and ordered trees with newick utilities (Junier, Thomas, and Evgeny M. Zdobnov. "The Newick utilities: high-throughput phylogenetic tree processing in the UNIX shell." Bioinformatics 26, no. 13 (2010): 1669-1670.), rooted and ordered trees saved in trees folder:
	nw_reroot nematodes.tre TSPI | nw_order -> nematodes_RO.pep
	nw_reroot fungi.tre EGLA | nw_order -> fungi_RO.pep
	nw_reroot insects.tre PHUM | nw_order -> insects_RO.pep

3) Ran tree_2_ultrametric.R script for each of the *_RO.pep trees. Ultrametrized-timed trees saved as *_RO_ultrametric.pep in trees folder. The R 3.5.1 script uses ape (Paradis E. & Schliep K. 2018. ape 5.0: an environment for modern phylogenetics and
  evolutionary analyses in R. Bioinformatics 35: 526-528.) and phytools (Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other
  things). Methods Ecol. Evol. 3 217-223. doi:10.1111/j.2041-210X.2011.00169.x) to time-calibrate and ultrametrize the an input tree with the 'makeChronosCalib' and 'chronos' functions.

	Divergence time for fungi: 241 mys
	Divergence time for insects: 358 mys
	Divergence time for nematodes: 428 mys

4) Compute gene counts per og per species table using CAFE_OGs_2_geneCounts.py script.

5) Obtained lambda tree for each of the 3 ultrametric rooted trees, to use later in CAFE for the lambda parameter. Replaced each node or leaf instance with '1' and saved file in trees folder.
	lambda_tree_*.txt

6) Ran 1st CAFE 4.2.1 analysis to find lambda and mu:

	~/Software/CAFE-4.2.1/CAFE/release/cafe run_insects.sh
	~/Software/CAFE-4.2.1/CAFE/release/cafe run_insects.sh 
	~/Software/CAFE-4.2.1/CAFE/release/cafe run_insects.sh 

	The bash script is the following (here showing fungi only, but performed the same on insects and nematodes):
"
	#!shell
	date

	# load ultrametric tree
	tree -i trees/fungi_RO_ultrametric.pep

	# load gene counts, filtering out non-rooting families
	load -i gene_counts/cafe_fungi.txt -p 0.05 -t 8 -l cafe_runs/output_fungi.txt

	# generate log
	log cafe_runs/log_fungi.txt

	# search for lambda
	lambdamu -s -t (((((1,1)1,1)1,((((1,1)1,1)1,1)1,((1,1)1,1)1)1)1,((1,1)1,((1,1)1,1)1)1)1,1) 

	# generate report
	report cafe_runs/report_run_fungi_p005
"

7) Ran 2nd CAFE 4.2.1 analysis to obtain the expansions and contractions tree annotations with the lambda and mu found before:

	~/Software/CAFE-4.2.1/CAFE/release/cafe run_insects.sh
	~/Software/CAFE-4.2.1/CAFE/release/cafe run_insects.sh 
	~/Software/CAFE-4.2.1/CAFE/release/cafe run_insects.sh

	The bash script is the following (here showing fungi only, but performed the same on insects and nematodes):
"
	#!shell
	date

	# load ultrametric tree
	tree -i trees/fungi_RO_ultrametric.pep

	# load gene counts, filtering out non-rooting families
	load -i gene_counts/cafe_fungi.txt -p 0.05 -t 8 -l cafe_runs/output_fungi.txt

	# generate log
	log cafe_runs/log_fungi.txt

	# search for lambda
	# lambdamu -s -t (((((1,1)1,1)1,((((1,1)1,1)1,1)1,((1,1)1,1)1)1)1,((1,1)1,((1,1)1,1)1)1)1,1) 
	lambdamu -l 0.00402051652068 -m 0.00428096307926 -t (((((1,1)1,1)1,((((1,1)1,1)1,1)1,((1,1)1,1)1)1)1,((1,1)1,((1,1)1,1)1)1)1,1)

	# generate report
	report cafe_runs/report_run_fungi_p005
"
	
