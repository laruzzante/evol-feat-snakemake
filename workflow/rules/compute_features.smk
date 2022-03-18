
rule compute_average_copy_number:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/ACN.tsv'
    script:
        '../scripts/compute_features/compute_ACN.py'
    log:
        'log/compute_average_copy_number.log'

rule compute_copy_number_variation:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/CNV.tsv'
    script:
        '../scripts/compute_features/compute_CNV.py'
    log:
        'log/compute_copy_number_variation.log'

rule compute_universality:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        info = rules.process_orthology_table.output.info
    output:
        'output/computed_orthogroup_features/UNI.tsv'
    script:
        '../scripts/compute_features/compute_UNI.py'
    log:
        'log/compute_universality.log'


rule compute_duplicability:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/DUP.tsv'
    script:
        '../scripts/compute_features/compute_DUP.py'
    log:
        'log/compute_duplicability.log'


rule compute_age:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_branchlengths = rules.create_MRCA_branchlengths_table.output[0]
    output:
        'output/computed_orthogroup_features/AGE.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/compute_features/compute_AGE.py'
    log:
        'log/compute_age.log'


rule compute_relative_universality:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0]
    output:
        'output/computed_orthogroup_features/RUN.tsv'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/compute_features/compute_RUN.py'
    log:
        'log/compute_relative_universality.log'


rule compute_expansions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        copy_number_variation_table = rules.create_copy_number_variation_table.output.copy_number_variation_table
    output:
        'output/computed_orthogroup_features/EXP.tsv'
    script:
        '../scripts/compute_features/compute_EXP.py'
    log:
        'log/compute_expansions.log'


rule compute_stabilities:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        copy_number_variation_table = rules.create_copy_number_variation_table.output.copy_number_variation_table
    output:
        'output/computed_orthogroup_features/STA.tsv'
    script:
        '../scripts/compute_features/compute_STA.py'
    log:
        'log/compute_stabilities.log'


rule compute_contractions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        copy_number_variation_table = rules.create_copy_number_variation_table.output.copy_number_variation_table
    output:
        'output/computed_orthogroup_features/CON.tsv'
    script:
        '../scripts/compute_features/compute_CON.py'
    log:
        'log/compute_contractions.log'


rule compute_relative_expansions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0],
        expansions = rules.compute_expansions.output
    output:
        'output/computed_orthogroup_features/REX.tsv'
    script:
        '../scripts/compute_features/compute_REX.py'
    log:
        'log/compute_relative_expansions.log'


rule compute_relative_stabilities:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0],
        stabilities = rules.compute_stabilities.output
    output:
        'output/computed_orthogroup_features/RST.tsv'
    script:
        '../scripts/compute_features/compute_RST.py'
    log:
        'log/compute_relative_stabilities.log'


rule compute_relative_contractions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0],
        contractions = rules.compute_contractions.output
    output:
        'output/computed_orthogroup_features/RCO.tsv'
    script:
        '../scripts/compute_features/compute_RCO.py'
    log:
        'log/compute_relative_contractions.log'


rule compute_synteny:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        synteny_counts = rules.create_synteny_counts_table.output.synteny_counts
    output:
        'output/computed_orthogroup_features/SYN.tsv'
    script:
        '../scripts/compute_features/compute_SYN.py'
    log:
        'log/compute_synteny.log'


rule compute_relative_synteny:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        synteny_counts = rules.create_synteny_counts_table.output.synteny_counts,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0]
    output:
        'output/computed_orthogroup_features/RSY.tsv'
    script:
        '../scripts/compute_features/compute_RSY.py'
    log:
        'log/compute_relative_synteny.log'


rule compute_maximum_synteny:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        synteny_counts = rules.create_synteny_counts_table.output.synteny_counts
    output:
        'output/computed_orthogroup_features/MSY.tsv'
    script:
        '../scripts/compute_features/compute_MSY.py'
    log:
        'log/compute_maximum_synteny.log'
