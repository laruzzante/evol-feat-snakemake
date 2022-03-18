MAX_MEMORY = config["MAX_MEMORY"]
MAX_RUNTIME = config["MAX_RUNTIME"] # in seconds
MAX_THREADS = config["MAX_THREADS"]

rule dimensionality_reductions:
    input:
        features = rules.merge_orthogroup_features.output[0]
    output:
        pca_coordinates = 'output/cluster_analysis/pca_coordinates.tsv',
        tsne_coordinates = 'output/cluster_analysis/tsne_coordinates.tsv',
        umap_coordinates = 'output/cluster_analysis/umap_coordinates.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/dimensionality_reductions.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/cluster_features/dimensionality_reductions.R'


rule homoscedasticity:
    input:
        features = rules.merge_orthogroup_features.output[0]
    output:
        plot = 'output/cluster_analysis/homoscedasticity.pdf'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/homoscedasticity.log'
    conda:
        '../envs/normality.yaml'
    script:
        '../scripts/cluster_features/homoscedasticity.R'


rule normality:
    input:
        features = rules.merge_orthogroup_features.output[0]
    output:
        plot = 'output/cluster_analysis/normality.pdf'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/normality.log'
    conda:
        '../envs/normality.yaml'
    script:
        '../scripts/cluster_features/normality.R'


rule pca:
    input:
        features = rules.merge_orthogroup_features.output[0]
    output:
        plot = 'output/cluster_analysis/pca.pdf'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME  # 6 hours = 21600 seconds
    log:
        'log/pca.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/cluster_features/pca.R'


# Hierarchichal Clustering with Pearson correlation distance and Average method
# Not using dimension reductions because I want to see the features dendrograms
rule hierarchichal_clustering:
    input:
        features = rules.merge_orthogroup_features.output[0]
    output:
        heatmap = 'output/cluster_analysis/hierarchichal_clustering.pdf'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/hierarchichal_clustering.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/cluster_features/hierarchichal_clustering.R'


# OPTICS clustering on tSNE
rule optics:
    input:
        tsne = rules.dimensionality_reductions.output.tsne_coordinates,
        features = rules.merge_orthogroup_features.output[0]
    output:
        plot = 'output/cluster_analysis/optics.pdf', # Plotting the optics memberships over umap coordinates
        optics_clusters = 'output/cluster_analysis/optics.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/optics.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/cluster_features/optics.R'


# DBSCAN clustering on variance-weighted PC
rule dbscan:
    input:
        umap = rules.dimensionality_reductions.output.umap_coordinates,
        features = rules.merge_orthogroup_features.output[0]
    output:
        weighted_pc_plot = 'output/cluster_analysis/weighted_pc_dbscan.pdf', # Plotting the dbscan memberships over tsne coordinates
        weighted_dbscan_clusters = 'output/cluster_analysis/weighted_pc_dbscan.tsv',
        umap_plot = 'output/cluster_analysis/umap_dbscan.pdf', # Plotting the dbscan memberships over tsne coordinates
        umap_dbscan_clusters = 'output/cluster_analysis/umap_dbscan.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/dbscan.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/cluster_features/dbscan.R'


# HDBSCAN clustering on UMAP
# rule hdbscan:
#     input:
#         umap = rules.dimensionality_reductions.output.umap_coordinates,
#         features = rules.merge_orthogroup_features.output[0]
#     output:
#         plot = 'output/cluster_analysis/hdbscan.pdf', # Plotting the hdbscan memberships over tsne coordinates
#         hdbscan_clusters = 'output/cluster_analysis/hdbscan.tsv'
#     threads: MAX_THREADS
#     resources:
#         mem_mb = MAX_MEMORY,  # 6 hours = 21600 seconds
#     log:
#         log = 'log/hdbscan.log'
#     conda:
#         '../envs/cluster_analysis.yaml'
#     script:
#         '../scripts/cluster_features/hdbscan.R'


# Not using dimension reductions because I want to see the features contributions
rule self_organising_map:
    input:
        features = rules.merge_orthogroup_features.output[0]
    output:
        plot = 'output/cluster_analysis/som.pdf', # SOM plot with R Kohonen map and counts_per_cell heatmap
        # som_clusters = 'output/cluster_analysis/som_clusters.tsv',
        # kmeans_som_superclusters = 'output/cluster_analysis/som_kmeans_superclusters.tsv', # superclusters based on kmeans
        # dbscan_som_superclusters = 'output/cluster_analysis/som_dbscan_superclusters.tsv' # superclusters based on DBSCAN
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/self_organising_map.log'
    conda:
        '../envs/self_organising_map.yaml'
    script:
        '../scripts/cluster_features/self_organising_map.R'


rule cluster_features:
    input:
        rules.dimensionality_reductions.output,
        rules.pca.output,
        # rules.hierarchichal_clustering.output,
        rules.optics.output,
        rules.dbscan.output,
        # rules.hdbscan.output,
        rules.self_organising_map.output
    output:
        touch('output/cluster_analysis/.done')
