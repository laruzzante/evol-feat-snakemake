rule format_OrthoDB_table:
    input:
        input_list['orthology_table']
    output:
        'output/formatted_orthology_table.tsv'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/format_OrthoDB_table.py'

# Might need to implement a Checkpoint instead, in order to get either OMA table or OrthoDB
# formatted table, depending on the format of the input orthoology table.
# For now, only OrthoDB formatting implemented.

rule process_orthology_table:
    input:
        rules.format_OrthoDB_table.output[0]
    output:
        orthogroups = 'output/orthogroups.pickle',
        genes = 'output/genes.pickle',
        species = 'output/species.pickle',
        info = 'output/orthology_info.txt'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/process_orthology_table.py'
