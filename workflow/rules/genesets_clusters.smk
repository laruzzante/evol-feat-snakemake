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

rule pvclust_on_sets:
    input:
        rules.get_orthogroup_sets_features.output[0]
    output:
        'output/genesets_cluster_analysis/hierarchical_clusters_of_orthogroup_sets.pdf'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/pvclust.R'

rule pca_on_sets:
    input:
        rules.get_orthogroup_sets_features.output[0]
    output:
        'output/genesets_cluster_analysis/pca_of_orthogroup_sets.pdf'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/pca.R'

rule pairwise_comparisons_on_sets:
    input:
        rules.get_orthogroup_sets_features.output[0]
    output:
        'output/genesets_cluster_analysis/pairwise_comparisons_of_orthogroup_sets.pdf'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/pairwise_comparisons.R'

rule tsne_on_sets:
    input:
        rules.get_orthogroup_sets_features.output[0]
    output:
        'output/genesets_cluster_analysis/tsne_of_orthogroup_sets.pdf'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/tsne.R'
