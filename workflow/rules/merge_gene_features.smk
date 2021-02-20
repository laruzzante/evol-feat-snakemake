
def get_user_gene_feature_files(wildcards):
    '''
    CODE COMMENTS
    '''
    user_gene_features_dir = checkpoints.split_user_gene_features.get().output[0] # get() here forces the checkpoint to rerun the DAG. E.g. without get(), I would only get a string of the output name.
    input_list = expand(os.path.join(user_gene_features_dir, '{feature}.tsv'),
                        feature=glob_wildcards(os.path.join(user_gene_features_dir, '{feature}.tsv')).feature)
    return input_list

rule merge_gene_features:
    input:
        get_user_gene_feature_files,
        expand('output/computed_gene_features/{feature}.tsv', feature=GENE_FEATURES_TO_COMPUTE)
    output:
        'output/merged_gene_features.tsv'
    script:
        '../scripts/merge_gene_features.py'
