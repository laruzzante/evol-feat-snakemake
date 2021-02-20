
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
        config["ultrametric_species_tree"]
    output:
        'output/mrca_branchlengths.tsv'
    script:
        '../scripts/create_mrca_branchlengths_table.R'


rule compute_age:
    input:
        rules.process_orthology_table.output.orthogroups,
        mrca_branchlengths = rules.create_mrca_branchlengths_table.output[0]
    output:
        'output/computed_orthogroup_features/AGE.tsv'
    script:
        '../scripts/compute_features/compute_AGE.py'
