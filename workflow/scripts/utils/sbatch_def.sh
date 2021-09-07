#!/bin/bash -l

#SBATCH --mail-type FAIL
#SBATCH --mail-user livio.ruzzante@unil.ch

#SBATCH --chdir /work/FAC/FBM/DEE/rwaterho/evofun/livio/software/CAFE5/bin
#SBATCH --job-name cafe5_orthogroup
#SBATCH --error=/work/FAC/FBM/DEE/rwaterho/evofun/livio/cafe_runs/badogsearch/err/cafe5_orthogroup.err
#SBATCH --output=/work/FAC/FBM/DEE/rwaterho/evofun/livio/cafe_runs/badogsearch/out/cafe5_orthogroup.out

#SBATCH --time=0-00:59:59

#SBATCH --cpus-per-task 48 
#SBATCH --mem 24G

source ~/.bashrc

conda activate cafe5

./cafe5 -i /work/FAC/FBM/DEE/rwaterho/evofun/livio/sbatch_scripts/cafe_badog_search/gene_count_tables/orthogroup_arthropoda_gene_counts.tsv -c 48 -t /work/FAC/FBM/DEE/rwaterho/evofun/livio/evol-feat-snakemake/workflow/input/arthropoda_rooted_time.tree -o /work/FAC/FBM/DEE/rwaterho/evofun/livio/cafe_runs/badogsearch/orthogroup > /work/FAC/FBM/DEE/rwaterho/evofun/livio/cafe_runs/badogsearch/log/orthogroup.log

