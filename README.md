# evol-feat-snakemake
Snakemake implementation of the evol-feat pipeline.

## Installation and usage
Clone or manually download the git repository:
git clone git@gitlab.com:evogenlab/evol-feat-snakemake.git

Go into the evol-feat-directory:
cd evol-feat-snakemake

Create conda environment from yaml environment file:
conda env create -f environment.yaml

Activate environment:
conda activate evol-feat

Go into the workflow directory:
cd workflow

Run snakemake to execute the workflow (either specify number of preferred cores or stick to 1):
snakemake --cores 1

Produce the DAG diagram:
snakemake --cores 1 --dag | dot -Tsvg > dag.svg



## Evolutionary Features
* ### UNI
universality: orthogroup species-span
* ### DUP
duplicability: proportion of species with gene duplicates
* ### ACN
average copy-number: average of per species gene copy-number
* ### CNV
copy-number variation: per species standard deviation of gene copies, divided by average copy-number