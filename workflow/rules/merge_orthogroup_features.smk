
def get_user_orthogroup_feature_files(wildcards):
    '''
    CODE COMMENTS
    '''
    user_orthogroup_features_dir = checkpoints.split_user_orthogroup_features.get().output[0] # get() here forces the checkpoint to rerun the DAG. E.g. without get(), I would only get a string of the output name.
    # user_orthogroup_feature_dir = checkpoints.split_user_orthogroup_features. # get() here forces the checkpoint to rerun the DAG. E.g. without get(), I would only get a string of the output name.
    # print(user_orthogroup_features_dir)
    input_list = expand(os.path.join(user_orthogroup_features_dir, '{feature}.tsv'),
                        feature=glob_wildcards(os.path.join(user_orthogroup_features_dir, '{feature}.tsv')).feature)
    # print(input_list)
    return input_list

rule merge_orthogroup_features:
    input:
        get_user_orthogroup_feature_files,
        expand('output/computed_orthogroup_features/{feature}.tsv', feature=OG_FEATURES_TO_COMPUTE)
        # get_og_feature_files(),
        # orthogroup_features = 'output/orthogroup_features/{feature}.tsv',
        # gene_features = 'output/gene_features/{feature}.tsv'
    output:
        'output/merged_orthogroup_features.tsv',
        # merged_gene_features = 'output/merged_gene_features.tsv'
    script:
        '../scripts/merge_orthogroup_features.py'
