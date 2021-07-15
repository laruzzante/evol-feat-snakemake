def get_input():

    input_list = {}
    input_list['user_orthogroup_features'] = []
    input_list['user_gene_features'] = []
    input_list['orthogroup_features_to_compute'] = []
    input_list['gene_features_to_compute'] = []
    input_list['all_orthogroup_features'] = []
    input_list['all_gene_features'] = []
    input_list['user_orthogroup_features_files'] = []
    input_list['user_gene_features_files'] = []
    input_list['orthology_table'] = ''
    input_list['ultrametric_species_tree'] = ''
    input_list['gff_file'] = ''
    input_list['focus_species_list'] = []
    input_list['focus_gene_list_files'] = []
    input_list['gene_sets'] = []
    USER_ORTHOGROUP_FEATURES = []
    USER_GENE_FEATURES = []

    if config['user_orthogroup_features_files']:
        for input_file in config['user_orthogroup_features_files']:
            with open(input_file) as infile:
                USER_ORTHOGROUP_FEATURES = infile.readline()[:-1].split('\t')[1:]
            input_list['user_orthogroup_features_files'].append(input_file)
            input_list['user_orthogroup_features'].append(USER_ORTHOGROUP_FEATURES)

    if config['user_gene_features_files']:
        for input_file in config['user_gene_features_files']:
            with open(input_file) as infile:
                USER_GENE_FEATURES = infile.readline()[:-1].split('\t')[1:]
            input_list['user_gene_features_files'].append(input_file)
            input_list['user_gene_features'].append(USER_GENE_FEATURES)

    if config['default_orthogroup_features']:
        SELECTED_ORTHOGROUP_FEATURES = config['selected_orthogroup_features']
        ORTHOGROUP_FEATURES_TO_COMPUTE = [feature for feature in SELECTED_ORTHOGROUP_FEATURES if feature not in USER_ORTHOGROUP_FEATURES]
        # print(f'Will compute following features: {ORTHOGROUP_FEATURES_TO_COMPUTE}')
        ORTHOGROUP_FEATURES_TO_NOT_COMPUTE = [feature for feature in SELECTED_ORTHOGROUP_FEATURES if feature in USER_ORTHOGROUP_FEATURES]
        # if ORTHOGROUP_FEATURES_TO_NOT_COMPUTE:
            # print(f'Warning: selected orthogroup features already provided by user, will not compute following: {ORTHOGROUP_FEATURES_TO_NOT_COMPUTE}')
        input_list['orthogroup_features_to_compute'] = ORTHOGROUP_FEATURES_TO_COMPUTE
    # else:
        # print(f'Warning: no DEFAULT ORTHOGROUP FEATURES to compute specified in the CONFIG file.')

    if config['default_gene_features']:
        SELECTED_GENE_FEATURES = config['selected_gene_features']
        GENE_FEATURES_TO_COMPUTE = [feature for feature in SELECTED_GENE_FEATURES if feature not in USER_GENE_FEATURES]
        # print(f'Will compute following features: {GENE_FEATURES_TO_COMPUTE}')
        GENE_FEATURES_TO_NOT_COMPUTE = [feature for feature in SELECTED_GENE_FEATURES if feature in USER_GENE_FEATURES]
        # if GENE_FEATURES_TO_NOT_COMPUTE:
            # print(f'Warning: selected gene features already provided by user, will not compute following: {GENE_FEATURES_TO_NOT_COMPUTE}')
        input_list['gene_features_to_compute'] = GENE_FEATURES_TO_COMPUTE
    # else:
    #     print(f'Warning: no DEFAULT GENE FEATURES to compute specified in the CONFIG file.')

    input_list['all_orthogroup_features'] = input_list['user_orthogroup_features'] + input_list['orthogroup_features_to_compute']
    input_list['all_gene_features'] = input_list['user_gene_features'] + input_list['gene_features_to_compute']

    if config['orthology_table']:
        input_list['orthology_table'] = config['orthology_table']

    if config['ultrametric_species_tree']:
        input_list['ultrametric_species_tree'] = config['ultrametric_species_tree']

    if config['gff_file']:
        input_list['gff_file'] = config['gff_file']

    if config['focus_species_list']:
        input_list['focus_species_list'] = config['focus_species_list']

    if config['focus_gene_list_files']:
        input_list['focus_gene_list_files'] = config['focus_gene_list_files']

    if config['gene_sets']:
        input_list['gene_sets'] = config['gene_sets']

    return input_list


# def get_missing_input_file_list(input_list):
#
#     missing_files = []
#
#     SET_1 = config['default_ortogroup_features']['set_1']
#     if any(input_list['orthogroup_features_to_compute']) == any(SET_1) and 'orthology_table' not in input_list.keys():
#         print(f'WARNING: missing ORTHOLOGY TABLE input file, cannot compute: {SET_1}')
#         missing_files.append('orthology_table')
#
#     SET_2 = config['default_ortogroup_features']['set_2']
#     if any(input_list['orthogroup_features_to_compute']) == any(SET_2) and ('orthology_table' not in input_list.keys() or 'ACN' not in input_list['all_orthogroup_features']):
#         if 'orthology_table' not in input_list.keys():
#             print(f'WARNING: missing ORTHOLOGY TABLE input file, cannot compute: {SET_2}')
#             missing_files.append('orthology_table')
#         if 'ACN' not in input_list['all_orthogroup_features']:
#             print(f'WARNING: missing ACN in user orthogroup features, cannot compute: {SET_2}')
#             missing_files.append('ACN feature')
#
#     SET_3 = config['default_ortogroup_features']['set_3']
#     if any(input_list['orthogroup_features_to_compute']) == any(SET_3) and 'ultrametric_species_tree' not in input_list.keys():
#         print(f'WARNING: missing UTLRAMETRIC SPECIES TREE input file, cannot compute: {SET_3}')
#         missing_files.append('ultrametric_species_tree')
#
#     return missing_files
