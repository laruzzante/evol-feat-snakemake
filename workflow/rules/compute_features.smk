# Here I tried to declare a "wildcarded rule name", without success:

# rule compute_feature:
#     input:
#         rules.create_orthology_dictionaries.output.orthogroups
#     output:
#         ex
#
#         expand("{dataset}/a.txt", dataset=DATASETS)
#         "aggregated.txt"
#     shell:
#         ...


rule compute_universality:
    input:
        rules.create_orthology_dictionaries.output.orthogroups
    output:
        UNI_orthogroups = 'results/evol_feat_orthogroups/UNI_orthogroups.tsv',
        UNI_genes = 'results/evol_feat_genes/UNI_genes.tsv'
    script:
        'compute_universality.py'

rule compute_duplicability:
    input:
        rules.create_orthology_dictionaries.output.orthogroups
    output:
        DUP_orthogroups = 'results/evol_feat_orthogroups/DUP_orthogroups.tsv',
        DUP_genes = 'results/evol_feat_genes/DUP_genes.tsv'
    script:
        'compute_duplicability.py'

rule compute_average_copy_number:
    input:
        rules.create_orthology_dictionaries.output.orthogroups
    output:
        ACN_orthogroups = 'results/evol_feat_orthogroups/ACN_orthogroups.tsv',
        ACN_genes = 'results/evol_feat_genes/ACN_genes.tsv'
    script:
        'compute_average_copy_number.py'

rule compute_copy_number_variation:
    input:
        rules.create_orthology_dictionaries.output.orthogroups,
        ACN_orthogroups = rules.compute_average_copy_number.output.ACN_orthogroups
    output:
        CNV_orthogroups = 'results/evol_feat_orthogroups/CNV_orthogroups.tsv',
        CNV_genes = 'results/evol_feat_genes/CNV_genes.tsv'
    script:
        'compute_copy_number_variation.py'

# rule compute_SYN:
#     input:
#         rules.create_orthology_dictionaries.output['orthogroups']
#         'data/{species}.gff'
#     output:
#         'results/{species}.tsv'
#     script:
#         'compute_UNI.py'
