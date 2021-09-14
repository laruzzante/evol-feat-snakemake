
rule compute_average_copy_number:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/ACN.tsv'
    script:
        '../scripts/compute_features/compute_ACN.py'

rule compute_copy_number_variation:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/CNV.tsv'
    script:
        '../scripts/compute_features/compute_CNV.py'


rule compute_universality:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        info = rules.process_orthology_table.output.info
    output:
        'output/computed_orthogroup_features/UNI.tsv'
    script:
        '../scripts/compute_features/compute_UNI.py'


rule compute_duplicability:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/DUP.tsv'
    script:
        '../scripts/compute_features/compute_DUP.py'


include: 'rules/preprocess_tree.smk'


rule compute_age:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        mrca_branchlengths = rules.create_mrca_branchlengths_table.output[0]
    output:
        'output/computed_orthogroup_features/AGE.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/compute_features/compute_AGE.py'


rule compute_relative_universality:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        mrca_ntips = rules.create_mrca_ntips_table.output[0]
    output:
        'output/computed_orthogroup_features/RUN.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/compute_features/compute_RUN.py'


include: 'rules/preprocess_cafe_results.smk'


rule compute_expansions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        copy_number_variation_table = rules.create_copy_number_variation_table.output.copy_number_variation_table
    output:
        'output/computed_orthogroup_features/EXP.tsv'
    script:
        '../scripts/compute_features/compute_EXP.py'


rule compute_relative_expansions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        copy_number_variation_table = rules.create_copy_number_variation_table.output.copy_number_variation_table,
        AGE = rules.compute_AGE.output[0]
    output:
        'output/computed_orthogroup_features/REX.tsv'
    script:
        '../scripts/compute_features/compute_REX.py'

# rule create_gene_counts_table:
#     input:
#         rules.process_orthology_table.output.orthogroups,
#         rules.process_orthology_table.output.species
#     output:
#         'output/gene_counts_table.tsv'
#     conda:
#         '../envs/basic.yaml'
#     script:
#         '../orthogroups_2_gene_counts.py'
#
#
# rule create_cafe_report:
#     input:
#         tree = input_list["ultrametric_species_tree"]
#         gene_counts = rules.create_gene_counts_table.output[0]
#     output:
#         'output/cafe_report.txt'
#     conda:
#         '../envs/cafe.yaml'
#     shell:
#         '
#         #!shell
#         date
#
#         # load ultrametric tree
#         tree -i trees/fungi_RO_ultrametric.pep
#
#         # load gene counts, filtering out non-rooting families
#         load -i gene_counts/cafe_fungi.txt -p 0.05 -t 8 -l cafe_runs/output_fungi.txt
#
#         # generate log
#         log cafe_runs/log_fungi.txt
#
#         # search for lambda
#         lambda -s
#
#         # correct for assembly copy-number errors
#         errormodel
#
#         # generate report
#         report cafe_runs/report_run_fungi_p005
#         '
