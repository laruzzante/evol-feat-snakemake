rule process_gff:
    input:
        genes = rules.process_orthology_table.output.genes,
        orthogroups = rules.process_orthology_table.output.orthogroups,
        gff = input_list["gff"]
    output:
        allspecies_synteny_scores = 'output/allspecies_synteny_scores.tsv'
    log:
        log = 'log/process_gff.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/process_gff.py'
