checkpoint split_user_orthogroup_features:
    input:
        input_list['user_orthogroup_features_files']
    output:
        directory('output/user_orthogroup_features/')
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/split_orthogroup_features.py'


checkpoint split_user_gene_features:
    input:
        input_list['user_gene_features_files']
    output:
        directory('output/user_gene_features/')
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/split_gene_features.py'
