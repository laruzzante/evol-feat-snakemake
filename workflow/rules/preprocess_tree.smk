rule create_mrca_branchlengths_table:
    input:
        input_list["ultrametric_species_tree"]
    output:
        'output/mrca_branchlengths.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/create_mrca_branchlengths_table.R'


rule create_mrca_ntips_table:
    input:
        input_list["ultrametric_species_tree"]
    output:
        'output/mrca_ntips.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/create_mrca_ntips_table.R'
