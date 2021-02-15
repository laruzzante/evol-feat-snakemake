rule process_orthology_table:
    input:
        config['orthology_table']
    output:
        orthogroups = 'output/orthogroups.pickle',
        genes = 'output/genes.pickle',
        info = 'output/orthology_info.txt'
    script:
        '../scripts/process_orthology_table.py'

rule compute_average_copy_number:
    input:
        rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/ACN.tsv'
        # ACN_genes = 'output/computed_gene_features/ACN_genes.tsv'
    script:
        '../scripts/compute_features/compute_ACN.py'

rule compute_copy_number_variation:
    input:
        rules.process_orthology_table.output.orthogroups,
        ACN_orthogroups = rules.compute_average_copy_number.output
    output:
        'output/computed_orthogroup_features/CNV.tsv'
        # CNV_genes = 'output/computed_gene_features/CNV_genes.tsv'
    script:
        '../scripts/compute_features/compute_CNV.py'


rule compute_universality:
    input:
        rules.process_orthology_table.output.orthogroups,
        info = rules.process_orthology_table.output.info
    output:
        'output/computed_orthogroup_features/UNI.tsv'
        # UNI_genes = 'output/computed_gene_features/UNI_genes.tsv'
    script:
        '../scripts/compute_features/compute_UNI.py'


rule compute_duplicability:
    input:
        rules.process_orthology_table.output.orthogroups
    output:
        'output/computed_orthogroup_features/DUP.tsv'
        # DUP_genes = 'output/computed_gene_features/DUP_genes.tsv'
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
        mrca_branchlengths = rules.create_mrca_branchlengths_table.output
    output:
        'output/computed_orthogroup_features/AGE.tsv'
        # AGE_genes = 'output/computed_gene_features/AGE_genes.tsv'
    script:
        '../scripts/compute_features/compute_AGE.py'
