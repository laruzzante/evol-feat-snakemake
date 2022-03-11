rule dimensionality_reductions:
    input:
        features = rules.merge_orthogroup_features.output[0]
    output:
        pca_coordinates = 'output/cluster_analysis/pca_coordinates.tsv',
        tsne_coordinates = 'output/cluster_analysis/tsne_coordinates.tsv',
        umap_coordinates = 'output/cluster_analysis/umap_coordinates.tsv'
    threads: 8
    resources:
        mem_mb = 128000,
        runtime_s = 21600  # 6 hours = 21600 seconds
    log:
        log = 'log/dimensionality_reductions.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/dimensionality_reductions.R'


# Hierarchichal Clustering with Pearson correlation distance and Average method
# Not using dimension reductions because I want to see the features dendrograms
rule hierarchichal_clustering:
    input:
        features = rules.merge_orthogroup_features.output[0]
    output:
        heatmap = 'output/cluster_analysis/hierarchichal_clustering.pdf'
    threads: 8
    resources:
        mem_mb = 128000,
        runtime_s = 21600  # 6 hours = 21600 seconds
    log:
        log = 'log/hierarchichal_clustering.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/hierarchichal_clustering.R'


# DBSCAN clustering on variance-weighted PC coordinates and on tSNE, both plotted on tSNE map
rule dbscan:
    input:
        pca_coordinates = rules.dimensionality_reductions.pca_coordinates
        tsne_coordinates = rules.dimensionality_reductions.tsne_coordinates
    output:
        plot = 'output/cluster_analysis/dbscan.pdf' # Simply plotting the dbscan memberships over tsne coordinates
        dbscan_clusters = 'output/cluster_analysis/dbscan.tsv'
        weighted_dbscan_clusters = 'output/cluster_analysis/weighted_dbscan.tsv'
    threads: 8
    resources:
        mem_mb = 128000,
        runtime_s = 21600  # 6 hours = 21600 seconds
    log:
        log = 'log/dbscan.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/dbscan.R'


# HDBSCAN clustering on on tSNE
rule hdbscan:
    input:
        tsne_coordinates = rules.dimensionality_reductions.tsne_coordinates
    output:
        plot = 'output/cluster_analysis/hdbscan.pdf' # Simply plotting the dbscan memberships over tsne coordinates
        hdbscan_clusters = 'output/cluster_analysis/hdbscan.tsv'
    threads: 8
    resources:
        mem_mb = 128000,
        runtime_s = 21600  # 6 hours = 21600 seconds
    log:
        log = 'log/hdbscan.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/hdbscan.R'


# OPTICS clustering on UMAP
rule optics:
    input:
        features = rules.dimensionality_reductions.umap_coordinates
    output:
        plot = 'output/cluster_analysis/optics.pdf' # Plotting the optics memberships over umap coordinates
        optics_clusters = 'output/cluster_analysis/optics.tsv'
    threads: 8
    resources:
        mem_mb = 128000,
        runtime_s = 21600  # 6 hours = 21600 seconds
    log:
        log = 'log/optics.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/optics.R'


# Not using dimension reductions because I want to see the features contributions
rule self_organising_map:
    input:
        features = rules.merge_orthogroup_features.output[0]
    output:
        plot = 'output/cluster_analysis/som.pdf', # SOM plot with R Kohonen map and counts_per_cell heatmap
        som_clusters = 'output/cluster_analysis/som_clusters.tsv',
        kmeans_som_superclusters = 'output/cluster_analysis/som_kmeans_superclusters.tsv', # superclusters based on kmeans
        dbscan_som_superclusters = 'output/cluster_analysis/som_dbscan_superclusters.tsv' # superclusters based on DBSCAN
    threads: 8
    resources:
        mem_mb = 128000,
        runtime_s = 21600  # 6 hours = 21600 seconds
    log:
        log = 'log/self_organising_map.log'
    conda:
        '../envs/cluster_analysis.yaml'
    script:
        '../scripts/self_organising_map.R'
