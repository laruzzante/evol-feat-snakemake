
rule extract_gene_lists_features:
    input:
        input_list['focus_gene_list_files'],
        rules.merge_gene_features.output[0]
    output:
        'output/gene_lists_features.tsv'
    log:
        'log/extract_gene_lists_features.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/extract_gene_lists_features.py'


rule get_orthogroup_features_by_gene:
    input:
        merged_orthogroup_features = rules.merge_orthogroup_features.output,
        orthogroups = rules.process_orthology_table.output.orthogroups,
        genes = rules.process_orthology_table.output.genes
    output:
        'output/orthogroup_features_by_gene.tsv'
    log:
        'log/get_orthogroup_features_by_gene.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/get_orthogroup_features_by_gene.py'


rule extract_focus_species_features:
    input:
        merged_orthogroup_features = rules.merge_orthogroups_features.output,
        orthogroup_features_by_gene = rules.get_orthogroup_features_by_gene,
        species = rules.process_orthology_table.output.species
    params:
        spec=FOCUS_SPECIES
    output:
        spec_orthogroups = expand('output/{spec}/merged_orthogroups_features.tsv', spec=FOCUS_SPECIES),
        spec_genes = expand('output/{spec}/orthogroup_features_by_gene.tsv', spec=FOCUS_SPECIES)
    log:
        'log/{spec}_extract_focus_species_features.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/extract_focus_species_features.py'
