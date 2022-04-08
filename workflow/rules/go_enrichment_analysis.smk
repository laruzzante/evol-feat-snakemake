MAX_MEMORY = config["MAX_MEMORY"]
MAX_RUNTIME = config["MAX_RUNTIME"] # in seconds
MAX_THREADS = config["MAX_THREADS"]


rule topgo_enrichment_on_som_clusters:
    input:
        som_clusters = rules.self_organising_map.output.som_clusters
        go_universe = input_list["go_universe"]
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
        'scripts/go_enrichment_analysis/topgo_enrichment_on_som_clusters.R'


rule topgo_enrichment_on_som_hc_superclusters:
    input:
        som_hc_superclusters = rules.self_organising_map.output.som_hc_superclusters
        go_universe = input_list["go_universe"]
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
        'scripts/go_enrichment_analysis/topgo_enrichment_on_som_hc_superclusters.R'
