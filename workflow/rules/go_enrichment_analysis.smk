MAX_MEMORY = config["MAX_MEMORY"]
MAX_RUNTIME = config["MAX_RUNTIME"] # in seconds
MAX_THREADS = config["MAX_THREADS"]


rule genes_go_2_orthogroups_go:
    input:
        genes = rules.process_orthology_table.output.genes,
        genes_go_universe = input_list["go_universe"]
    output:
        orthogroups_go_universe = 'output/go_enrichment_analysis/orthogroups_go_universe.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/genes_go_2_orthogroups_go.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/go_enrichment_analysis/genes_go_2_orthogroups_go.py'


rule topgo_enrichment_on_feature_ranks:
    input:
        features = rules.merge_orthogroup_features.output,
        go_universe = rules.genes_go_2_orthogroups_go.output.orthogroups_go_universe
    output:
        feature_ranks_go = 'output/go_enrichment_analysis/feature_ranks_go.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/topgo_enrichment_on_feature_ranks.log'
    conda:
        '../envs/go_enrichment_analysis.yaml'
    script:
        '../scripts/go_enrichment_analysis/topgo_enrichment_on_feature_ranks.R'


rule topgo_enrichment_on_som_clusters:
    input:
        som_clusters = rules.self_organising_map.output.som_clusters,
        go_universe = rules.genes_go_2_orthogroups_go.output.orthogroups_go_universe
    output:
        som_clusters_go = 'output/go_enrichment_analysis/som_clusters_go.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/topgo_enrichment_on_som_clusters.log'
    conda:
        '../envs/go_enrichment_analysis.yaml'
    script:
        '../scripts/go_enrichment_analysis/topgo_enrichment_on_som_clusters.R'


rule topgo_enrichment_on_som_hc_superclusters:
    input:
        som_hc_superclusters = rules.self_organising_map.output.som_hc_superclusters,
        go_universe = rules.genes_go_2_orthogroups_go.output.orthogroups_go_universe
    output:
        som_hc_superclusters_go = 'output/go_enrichment_analysis/som_hc_superclusters_go.tsv'
    threads: MAX_THREADS
    resources:
        mem_mb = MAX_MEMORY,
        runtime_s = MAX_RUNTIME
    log:
        'log/topgo_enrichment_on_som_hc_superclusters.log'
    conda:
        '../envs/go_enrichment_analysis.yaml'
    script:
        '../scripts/go_enrichment_analysis/topgo_enrichment_on_som_hc_superclusters.R'
