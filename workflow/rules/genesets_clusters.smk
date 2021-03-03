rule genesets_2_ogsets:
    input:
        input_list['gene_sets'],
        genes = rules.process_orthology_table.output.genes
    output:
        'output/orthogroup_sets.tsv'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/genesets_2_orthogroupsets.py'

rule get_orthogroup_sets_features:
    input:
        ogsets = rules.genesets_2_ogsets.output[0],
        features = rules.merge_orthogroup_features.output[0]
    output:
        'output/orthogroup_sets_features.tsv'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/get_orthogroup_sets_features.py'

rule hclust_sets:
    input:
        rules.get_orthogroup_sets_features.output[0]
    output:
        'output/cluster_analysis/hierarchical_clusters_orthogroup_sets.pdf'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/cluster_analysis/hclust_sets.R'
