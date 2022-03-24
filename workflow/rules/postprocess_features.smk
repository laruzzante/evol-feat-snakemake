# def get output filenames from input focus species


# rule extract_focus_species_features:
#     input:
#         rules.merge_gene_features.output[0]
#     params:
#         spec=FOCUS_SPECIES
#     output:
#         expand('output/{spec}_gene_features.tsv', spec=FOCUS_SPECIES)
#     log:
#         'log/extract_focus_species_features.log'
#     conda:
#         '../envs/basic.yaml'
#     run:
#         shell("grep {params.spec}: {input}>>{output}")


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
        orthogroups = rules.process_orthology_table.orthogroups,
        genes = rules.process_orthology_table.genes
    output:
        'output/orthogroup_features_by_gene.tsv'
    log:
        'log/get_orthogroup_features_by_gene.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/get_orthogroup_features_by_gene.py'
