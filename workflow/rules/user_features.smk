checkpoint split_user_orthogroup_features:
    input:
        config['user_orthogroup_features_files']
    output:
        directory('output/user_orthogroup_features/')
    script:
        '../scripts/split_orthogroup_features.py'
