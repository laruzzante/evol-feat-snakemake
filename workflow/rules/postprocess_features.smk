# def get output filenames from input focus species


rule extract_focus_species_features:
    input:
        rules.merge_gene_features.output[0]
    params:
        spec=FOCUS_SPECIES
    output:
        expand('output/{spec}_gene_features.tsv', spec=FOCUS_SPECIES)
    run:
        shell("grep {params.spec}: {input}>>{output}")


rule extract_gene_lists_features:
    input:
        input_list['focus_gene_list_files'],
        rules.merge_gene_features.output[0]
    output:
        'output/gene_lists_features.tsv'
    script:
        '../scripts/extract_gene_lists_features.py'
