# evol-feat-snakemake
Snakemake implementation of the evol-feat pipeline.

## Installation and usage
* Clone or manually download the git repository:
```bash
git clone git@gitlab.com:evogenlab/evol-feat-snakemake.git
```

* Go into the evol-feat-directory:
```bash
cd evol-feat-snakemake
```

* Create conda environment from yaml environment file:
```bash
conda env create -f environment.yaml
```

* Activate environment:
```bash
conda activate evol-feat
```

* Go into the workflow directory:
```bash
cd workflow
```

* Run snakemake to execute the workflow (either specify the number of preferred cores or stick to 1):
```bash
snakemake --cores 1
```

* Produce the Directed Acyclic Graph (DAG) diagram:
```bash
snakemake --cores 1 --dag | dot -Tsvg > dag.svg
```

## Evolutionary Features
* ### UNI
universality: orthogroup species-span
* ### DUP
duplicability: proportion of species with gene duplicates
* ### ACN
average copy-number: average of per species gene copy-number
* ### CNV
copy-number variation: per species standard deviation of gene copies, divided by average copy-number
