rule create_copy_number_variation_table:
    input:
        cafe_results = input_list['cafe_results']
    output:
        copy_number_variation_table = 'output/copy_number_variation_table.tsv'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/process_cafe5_base_change.py'
