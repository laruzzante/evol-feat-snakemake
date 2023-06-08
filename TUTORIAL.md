Evol-Feat tutorial:

## Orthology Table

1. Download the OG2genes tab flat file from https://data.orthodb.org/download/

2. **Orthogroups to genes** relationships are defined in OrthoDB by taxonomic level hierarchies. Filter out the genes at the taxonomic level of interest with the following bash command (example given for the Apoidea (bees) level: NCBI taxid = 34735), and save the output into a text file (after the '>' pipe). Replace 'at34735\t' with 'atXXXX\t' and the taxonomic level of preference:
```bash
	zcat odb11v0_OG2genes.tab.gz | sed -n 'at34735\t' > apoidea_odb10v1_OG2genes.tab
```

## Species Tree

3. Have your **ultrametric species tree** ready (species names must correspond to the species ids used in the odb_OG2gene file). Otherwise create it from OrthoDB using Orthophile:
	a. Extract all unique species ids at the taxonomic level of interest from the apoidea_odb10v1_OG2gene.tab file with the 'get_species_ids_from_odb_OG2genes.py' python script from the "scripts/utils" folder:
```bash
		python3 get_species_ids_from_odbOG2genes.py apoidea_odb11v0_OG2genes.tab
```
	b. Save the species id list into a text file (one id per row), and use it as input file for the Orthophile phylogenetic tree building workflow.
	c. Orthophile also works with NCBItaxid alone, but its easier to read a tree with the corresponding species names. You can thus edit your species list file by adding the species names before the taxid, separated by a comma. E.g. "103933_0, Bombus_bifarius". It is important to have underscores '_' instead of spaces in the species names between genus and species. If you instead have too many species and want to automate this step, run the script in "script/utils" named "get_NCBI_taxname.py". Note that this script accesses the NCBI Taxonomy database, and any "_X" after the taxid as formatted by OrthoDB must be removed. The command runs as follow:
```bash
		python3 get_NCBI_taxname.py taxid1 taxid2 taxid3 ...
```
