rule create_cafe_base_change_table:
    input:
        input_list['cafe_results']
    output:
        'output/cafe_base_change_table.tsv'
    conda:
        '../envs/basic.yaml'
    script:
        'process_cafe5_base_change.py'
