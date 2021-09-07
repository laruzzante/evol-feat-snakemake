
rule compute_average_copy_number:
    input:
        rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/ACN.tsv'
    script:
        '../scripts/compute_features/compute_ACN.py'

rule compute_copy_number_variation:
    input:
        rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/CNV.tsv'
    script:
        '../scripts/compute_features/compute_CNV.py'


rule compute_universality:
    input:
        rules.process_orthology_table.output.orthogroups,
        info = rules.process_orthology_table.output.info
    output:
        'output/computed_orthogroup_features/UNI.tsv'
    script:
        '../scripts/compute_features/compute_UNI.py'


rule compute_duplicability:
    input:
        rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/DUP.tsv'
    script:
        '../scripts/compute_features/compute_DUP.py'


rule create_mrca_branchlengths_table:
    input:
        input_list["ultrametric_species_tree"]
    output:
        'output/mrca_branchlengths.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/create_mrca_branchlengths_table.R'


rule compute_age:
    input:
        rules.process_orthology_table.output.orthogroups,
        mrca_branchlengths = rules.create_mrca_branchlengths_table.output[0]
    output:
        'output/computed_orthogroup_features/AGE.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/compute_features/compute_AGE.py'


rule create_mrca_ntips_table:
    input:
        input_list["ultrametric_species_tree"]
    output:
        'output/mrca_ntips.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/create_mrca_ntips_table.R'


rule compute_relative_universality:
    input:
        rules.process_orthology_table.output.orthogroups,
        mrca_ntips = rules.create_mrca_ntips_table.output[0]
    output:
        'output/computed_orthogroup_features/RUN.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/compute_features/compute_RUN.py'


rule compute_age:
    input:
        rules.process_orthology_table.output.orthogroups,
        mrca_branchlengths = rules.create_mrca_branchlengths_table.output[0]
    output:
        'output/computed_orthogroup_features/AGE.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/compute_features/compute_AGE.py'

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
