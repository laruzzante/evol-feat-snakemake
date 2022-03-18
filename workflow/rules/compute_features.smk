
rule compute_average_copy_number:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/ACN.tsv'
    log:
        'log/compute_average_copy_number.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_ACN.py'


rule compute_copy_number_variation:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/CNV.tsv'
    log:
        'log/compute_copy_number_variation.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_CNV.py'


rule compute_universality:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        info = rules.process_orthology_table.output.info
    output:
        'output/computed_orthogroup_features/UNI.tsv'
    log:
        'log/compute_universality.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_UNI.py'


rule compute_duplicability:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/DUP.tsv'
    log:
        'log/compute_duplicability.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_DUP.py'


rule compute_age:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_branchlengths = rules.create_MRCA_branchlengths_table.output[0]
    output:
        'output/computed_orthogroup_features/AGE.tsv'
    log:
        'log/compute_age.log'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/compute_features/compute_AGE.py'


rule compute_relative_universality:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0]
    output:
        'output/computed_orthogroup_features/RUN.tsv'
    log:
        'log/compute_relative_universality.log'
    conda:
        '../envs/phylogeny.yaml'
    script:
        '../scripts/compute_features/compute_RUN.py'


rule compute_expansions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        copy_number_variation_table = rules.create_copy_number_variation_table.output.copy_number_variation_table
    output:
        'output/computed_orthogroup_features/EXP.tsv'
    log:
        'log/compute_expansions.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_EXP.py'


rule compute_stabilities:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        copy_number_variation_table = rules.create_copy_number_variation_table.output.copy_number_variation_table
    output:
        'output/computed_orthogroup_features/STA.tsv'
    log:
        'log/compute_stabilities.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_STA.py'


rule compute_contractions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        copy_number_variation_table = rules.create_copy_number_variation_table.output.copy_number_variation_table
    output:
        'output/computed_orthogroup_features/CON.tsv'
    log:
        'log/compute_contractions.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_CON.py'


rule compute_relative_expansions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0],
        expansions = rules.compute_expansions.output
    output:
        'output/computed_orthogroup_features/REX.tsv'
    log:
        'log/compute_relative_expansions.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_REX.py'


rule compute_relative_stabilities:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0],
        stabilities = rules.compute_stabilities.output
    output:
        'output/computed_orthogroup_features/RST.tsv'
    log:
        'log/compute_relative_stabilities.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_RST.py'


rule compute_relative_contractions:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0],
        contractions = rules.compute_contractions.output
    output:
        'output/computed_orthogroup_features/RCO.tsv'
    log:
        'log/compute_relative_contractions.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_RCO.py'


rule compute_synteny:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        synteny_counts = rules.create_synteny_counts_table.output.synteny_counts
    output:
        'output/computed_orthogroup_features/SYN.tsv'
    log:
        'log/compute_synteny.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_SYN.py'


rule compute_relative_synteny:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        synteny_counts = rules.create_synteny_counts_table.output.synteny_counts,
        MRCA_ntips = rules.create_MRCA_ntips_table.output[0]
    output:
        'output/computed_orthogroup_features/RSY.tsv'
    log:
        'log/compute_relative_synteny.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_RSY.py'


rule compute_maximum_synteny:
    input:
        orthogroups = rules.process_orthology_table.output.orthogroups,
        synteny_counts = rules.create_synteny_counts_table.output.synteny_counts
    output:
        'output/computed_orthogroup_features/MSY.tsv'
    log:
        'log/compute_maximum_synteny.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/compute_features/compute_MSY.py'
