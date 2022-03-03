rule format_orthology_table:
    input:
        orthology_table = input_list['orthology_table']
    output:
        formatted_orthology_table = 'output/formatted_orthology_table.tsv'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/format_orthology_table.py'

# Might need to implement a Checkpoint instead, in order to get either OMA table or OrthoDB
# formatted table, depending on the format of the input orthology table.
# For now, only OrthoDB formatting implemented.

rule process_orthology_table:
    input:
        formatted_orthology_table = rules.format_orthology_table.output.formatted_orthology_table
    output:
        orthogroups = 'output/.orthogroups.pickle',
        genes = 'output/.genes.pickle',
        species = 'output/.species.pickle',
        orthogroups_2_species_2_genes = 'output/.orthogroups_2_species_2_genes.pickle',
        info = 'output/orthology_info.txt'
    log:
        log = 'log/process_orthology_table.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/process_orthology_table.py'
