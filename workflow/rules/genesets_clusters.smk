MAX_MEMORY = config["MAX_MEMORY"]
MAX_RUNTIME = config["MAX_RUNTIME"] # in seconds
MAX_THREADS = config["MAX_THREADS"]


rule genesets_2_ogsets:
    input:
        input_list['gene_sets'],
        genes = rules.process_orthology_table.output.genes
    output:
        'output/genesets_cluster_analysis/orthogroup_sets.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/genesets_2_ogsets.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/genesets_2_orthogroupsets.py'


rule get_orthogroup_sets_features:
    input:
        ogsets = rules.genesets_2_ogsets.output[0],
        features = rules.merge_orthogroup_features.output[0]
    output:
        'output/genesets_cluster_analysis/orthogroup_sets_features.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/get_orthogroup_sets_features.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/get_orthogroup_sets_features.py'


rule hierarchichal_clustering_on_sets:
    input:
        rules.get_orthogroup_sets_features.output[0]
    output:
        'output/genesets_cluster_analysis/hierarchical_clusters_of_orthogroup_sets.pdf'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/hierarchichal_clustering_on_sets.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/hierarchichal_clustering_on_sets.R'


rule pca_on_sets:
    input:
        rules.get_orthogroup_sets_features.output[0]
    output:
        'output/genesets_cluster_analysis/pca_of_orthogroup_sets.pdf'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/pca_on_sets.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/pca_on_sets.R'


rule pairwise_comparisons_on_sets:
    input:
        rules.get_orthogroup_sets_features.output[0]
    output:
        'output/genesets_cluster_analysis/pairwise_comparisons_of_orthogroup_sets.pdf'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/pairwise_comparisons_on_sets.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/pairwise_comparisons_on_sets.R'


rule dimensionality_reductions_on_sets:
    input:
        rules.get_orthogroup_sets_features.output[0]
    output:
        'output/genesets_cluster_analysis/dimensionality_reductions_of_orthogroup_sets.pdf'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/dimensionality_reductions_on_sets.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/genesets_cluster_analysis/dimensionality_reductions_on_sets.R'


rule genesets_clusters:
    input:
        rules.hierarchichal_clustering_on_sets.output,
        rules.pca_on_sets.output,
        rules.pairwise_comparisons_on_sets.output,
        rules.dimensionality_reductions_on_sets.output
    output:
        touch('output/genesets_cluster_analysis/.done')
