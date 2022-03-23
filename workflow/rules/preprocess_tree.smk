rule create_MRCA_branchlengths_table:
    input:
        input_list["ultrametric_species_tree"]
    output:
        'output/MRCA_branchlengths.tsv'
    log:
        'log/create_MRCA_branchlengths_table.log'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/create_MRCA_branchlengths_table.R'


rule create_MRCA_ntips_table:
    input:
        input_list["ultrametric_species_tree"]
    output:
        'output/MRCA_ntips.tsv'
    log:
        'log/create_MRCA_ntips_table.log'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/create_MRCA_ntips_table.R'
