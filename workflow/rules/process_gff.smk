MAX_MEMORY = config["MAX_MEMORY"]
MAX_RUNTIME = config["MAX_RUNTIME"] # in seconds
MAX_THREADS = config["MAX_THREADS"]

rule create_ordered_gff_genes_table:
    input:
        genes = rules.process_orthology_table.output.genes,
        orthogroups = rules.process_orthology_table.output.orthogroups,
        gff = input_list["gff"]
    output:
        ordered_gff_genes = 'output/ordered_gff_genes.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/create_ordered_gff_genes_table.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/create_ordered_gff_genes_table.py'


rule create_synteny_counts_table:
    input:
        orthogroups_2_species_2_genes = rules.process_orthology_table.output.orthogroups_2_species_2_genes,
        ordered_gff_genes = rules.create_ordered_gff_genes_table.output.ordered_gff_genes
    output:
        synteny_dict = 'output/.synteny.pickle',
        synteny_counts = 'output/synteny_counts.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/create_synteny_counts_table.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/create_synteny_counts_table.py'
