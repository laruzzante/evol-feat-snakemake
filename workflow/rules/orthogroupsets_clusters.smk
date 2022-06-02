
rule format_og_sets:
    input:
        input_list['orthogroup_sets'],
        orthogroups = rules.process_orthology_table.output.orthogroups
    output:
        'output/orthogroupsets_cluster_analysis/orthogroup_sets.tsv'
    log:
        'log/format_og_sets.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/format_orthogroupsets.py'


rule get_og_sets_features:
    input:
        ogsets = rules.format_og_sets.output[0],
        features = rules.merge_orthogroup_features.output[0]
    output:
        'output/orthogroupsets_cluster_analysis/orthogroup_sets_features.tsv'
    log:
        'log/get_og_sets_features.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/get_orthogroup_sets_features.py'


rule hierarchichal_clustering_on_og_sets:
    input:
        rules.get_og_sets_features.output[0]
    output:
        'output/orthogroupsets_cluster_analysis/hierarchical_clusters_of_orthogroup_sets.pdf'
    log:
        'log/hierarchichal_clustering_on_og_sets.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/hierarchichal_clustering_on_sets.R'


rule pca_on_og_sets:
    input:
        rules.get_og_sets_features.output[0]
    output:
        'output/orthogroupsets_cluster_analysis/pca_of_orthogroup_sets.pdf'
    log:
        'log/pca_on_og_sets.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/pca_on_sets.R'


rule pairwise_comparisons_on_og_sets:
    input:
        rules.get_og_sets_features.output[0]
    output:
        'output/orthogroupsets_cluster_analysis/pairwise_comparisons_of_orthogroup_sets.pdf'
    log:
        'log/pairwise_comparisons_on_og_sets.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/pairwise_comparisons_on_sets.R'


rule dimensionality_reductions_on_og_sets:
    input:
        rules.get_og_sets_features.output[0]
    output:
        'output/orthogroupsets_cluster_analysis/dimensionality_reductions_of_orthogroup_sets.pdf'
    log:
        'log/dimensionality_reductions_on_og_sets.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/dimensionality_reductions_on_sets.R'


rule orthogroupsets_clusters:
    input:
        rules.hierarchichal_clustering_on_og_sets.output,
        # rules.pca_on_og_sets.output,
        # rules.pairwise_comparisons_on_og_sets.output,
        # rules.dimensionality_reductions_on_og_sets.output
    output:
        touch('output/orthogroupsets_cluster_analysis/.done')
