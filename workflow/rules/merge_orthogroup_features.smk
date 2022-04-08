
def get_user_orthogroup_feature_files(wildcards):
    '''
    Function to parse full list of user provided orthogroup features files
    '''
    user_orthogroup_features_dir = checkpoints.split_user_orthogroup_features.get().output[0] # get() here forces the checkpoint to rerun the DAG. E.g. without get(), I would only get a string of the output name.
    input = expand(os.path.join(user_orthogroup_features_dir, '{feature}.tsv'),
                        feature=glob_wildcards(os.path.join(user_orthogroup_features_dir, '{feature}.tsv')).feature)
    return input

rule merge_orthogroup_features:
    input:
        get_user_orthogroup_feature_files,
        expand('output/computed_orthogroup_features/{feature}.tsv', feature=OG_FEATURES_TO_COMPUTE)
    output:
        'output/merged_orthogroup_features.tsv'
    log:
        'log/merge_orthogroup_features.log'
    conda:
        '../envs/basic.yaml'
    script:
        '../scripts/merge_orthogroup_features.py'
