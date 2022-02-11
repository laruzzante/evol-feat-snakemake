Pre-processing input files for SYNTENY (SYN) metric computation:

Formatting input files to get matching gene ids across orthology table (OrthoDB id) and GFFs ().

Download files from OrthoDB:

1) odb10v1_OG2genes.tab
2) odb10v1_genes.tab

Download GFF files

## Check which gene id to use, either VBid or ODBid

Extract orthogroups mapped at arthropoda node (i.e. taxid = 6656):
	grep -w 'at6656' odb10v1_OG2genes.tab > arthropoda_odb10v1_OG2genes.tab

We need to only consider geneids that are present in the arthropoda_odb10v1_OG2gene.tab, as the odb10v1_genes.tab
contains every gene in orthodb, not just arthropoda. We filter out these rows within the python script
after parsing and storing in a dictionary the arthropoda_odb10v1_OG2genes.tab ODB gene ids.

No mismatches were found between VBids and ODBids, i.e. each gene id had a 1:1 correspondance.

Hence, we can map the GFF VBids to ODBids and continue using ODBids along the pipeline, as it is the
default naming scheme used by the first input file, i.e. the OrthoDB orthology table.

I will use the 3rd column in odb10v1_genes.tab as VectorBase gene id, as this is the one corresponding to the ids in the GFFs.


## Converting ids in GFFs
Each GFF now has to be converted. But first we can speed up things by extracting CDS rows, with a grep command:

grep -w 'CDS' *for_every_GFF_file* > GFF_file_CDS

then  ran script gff_vbid_2_odbid.py to create the final converted GFF file with ODB gene ids.



Ran the CDS extraction and ID conversion for all GFFs with the 'grep_cds_and_convert_id.sh' script.

