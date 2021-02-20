
rule extract_focus_species_features:
    input:
        rules.merge_gene_features.output[0]
    params:
        spec=FOCUS_SPECIES_LIST
    output:
        expand('output/{spec}_gene_features.tsv', spec=FOCUS_SPECIES_LIST)
    run:
        # print(FOCUS_SPECIES_LIST)
        # shell("head -n 1 {input}>{output}")
        shell("grep {params.spec}: {input}>>{output}")
