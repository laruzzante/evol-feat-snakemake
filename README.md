Snakemake implementation of the evol-feat pipeline.

## Installation and usage
* Make sure to have the **conda** or **miniconda** environment manager system installed:
https://docs.conda.io/en/latest/miniconda.html

* Make sure to have git and snakemake activated:
https://git-scm.com/
https://snakemake.readthedocs.io/en/stable/

* Create and activate your conda snakemake environment:
```bash
conda create --name snakemake snakemake=7.3.2
conda activate snakemake
```

* Clone or manually download the git repository:
```bash
git clone git@gitlab.com:evogenlab/evol-feat-snakemake.git
```

* Move to the evol-feat-snakemake workflow directory:
```bash
cd evol-feat-snakemake/workflow
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
snakemake --dag | dot -Tsvg > dag.svg
```

## Computable Evolutionary Features
* ### UNI
UNIversality: The proportion of the total species present in an OG (all species, UNI=1) [requires an orthology table]
* ### DUP
DUPlicability: The proportion of species present in an OG that have multi-copy orthologues [requires an orthology table]
* ### ACN
Average Copy-Number: The average (mean) orthologue copy number across all species present in an OG [requires an orthology table]
* ### CNV
Copy-Number Variation: The standard deviation of orthologue counts per species present in an OG divided by the ACN [requires an orthology table]
* ### AGE
taxonomic AGE: Age of the last common ancestor of species in an OG, in terms of millions of years since divergence, computed from the ultrametric species phylogeny  [requires an orthology table and a time-calibrated phylogeny]
* ### RUN
Relative UNiversality: UNI divided by the number of nodes derived from OG’s last common ancestor [requires an orthology table and a time-calibrated phylogeny]
* ### EXP
EXPansions: CAFE quantified proportions of gene gain nodes for an OG [requires a CAFE summary report named base_asr from CAFE5]
* ### CON
CONtractions: CAFE quantified proportions of gene loss nodes for an OG [requires a CAFE summary report named base_asr from CAFE5]
* ### STA
STAbility: CAFE quantified proportions of no copy-number change nodes for an OG [requires a CAFE summary report named base_asr from CAFE5]
* ### SYN
SYNteny: The species-averaged proportion of orthologues in an OG that maintain their orthologous neighbours in the genomes of the other species [requires an orthology table and single-file concatenated per-species GFF files]
* ### RUN
Relative UNiversality: UNI divided by the number of nodes derived from OG’s last common ancestor [requires an orthology table and a time-calibrated phylogeny]
* ### REX
Relative EXpansions: EXP divided by the number of nodes derived from OG’s last common ancestor [requires a CAFE summary report named base_asr from CAFE5 a time-calibrated phylogeny]
* ### RCO
Relative COntractions: CON divided by the number of nodes derived from OG’s last common ancestor [requires a CAFE summary report named base_asr from CAFE5 a time-calibrated phylogeny]
* ### RST
Relative STability: STA divided by the number of nodes derived from OG’s last common ancestor [requires a CAFE summary report named base_asr from CAFE5 a time-calibrated phylogeny]
* ### RSY
Relative SYnteny: SYN divided by the number of nodes derived from OG’s last common ancestor [requires an orthology table, a single-file concatenated per-species GFF files and a time-calibrated phylogeny]
* ### MSY
Maximum SYnteny: The species-maximum proportion of orthologues in an OG that maintain their orthologous neighbours in the genomes of the other species  [requires an orthology table, a single-file concatenated per-species GFF files and a time-calibrated phylogeny]

## Suggested user-provided metrics:
* ### EVR
EVolutionary Rate: The average rate of protein sequence divergence normalised by the distance (% identity) between each pair of species as computed by OrthoDB [requires an OrthoDB Evolutionary Rates file]
