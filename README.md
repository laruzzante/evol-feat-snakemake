Snakemake implementation of the evol-feat pipeline.

## Installation and usage
* Make sure to have the **conda** or **miniconda** environment manager system installed:
https://docs.conda.io/en/latest/miniconda.html

* Clone or manually download the git repository:
```bash
git clone git@gitlab.com:evogenlab/evol-feat-snakemake.git
```

* Go into the evol-feat-directory:
```bash
cd evol-feat-snakemake
```

* Go into the workflow directory:
```bash
cd workflow
```

* Edit the config.yaml file:
```bash
vim config.yaml
```

* Rune the pipeline, where N is the number of cores you want snakemake to use:
```bash
snakemake --cores <N> --use-conda
```

* Produce the Directed Acyclic Graph (DAG) diagram:
```bash
snakemake --cores 1 --dag | dot -Tsvg > dag.svg
```

## Computable Evolutionary Features
* ### UNI
universality: orthogroup species-span
* ### DUP
duplicability: proportion of species with gene duplicates
* ### ACN
average copy-number: average of per species gene copy-number
* ### CNV
copy-number variation: per species standard deviation of gene copies, divided by average copy-number
* ### AGE
orthogroup age: age of the most recent common ancestor of the orthogroup's species span
* ### RUN
relative universality: orthogroup species-span relative to the clade originating from the orthogroup's last common ancestor
