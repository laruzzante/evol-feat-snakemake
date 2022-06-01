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
    input_list['cafe_results'] = ''
    input_list['gff'] = ''
    input_list['focus_species'] = []
    input_list['focus_gene_list_files'] = []
    input_list['gene_sets'] = []
    input_list['orthogroup_sets'] = []
    input_list['go_universe'] = []
    input_list['ontology'] = []
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
        ORTHOGROUP_FEATURES_TO_NOT_COMPUTE = [feature for feature in SELECTED_ORTHOGROUP_FEATURES if feature in USER_ORTHOGROUP_FEATURES]
        input_list['orthogroup_features_to_compute'] = ORTHOGROUP_FEATURES_TO_COMPUTE

    if config['default_gene_features']:
        SELECTED_GENE_FEATURES = config['selected_gene_features']
        GENE_FEATURES_TO_COMPUTE = [feature for feature in SELECTED_GENE_FEATURES if feature not in USER_GENE_FEATURES]
        GENE_FEATURES_TO_NOT_COMPUTE = [feature for feature in SELECTED_GENE_FEATURES if feature in USER_GENE_FEATURES]
        input_list['gene_features_to_compute'] = GENE_FEATURES_TO_COMPUTE

    input_list['all_orthogroup_features'] = input_list['user_orthogroup_features'] + input_list['orthogroup_features_to_compute']
    input_list['all_gene_features'] = input_list['user_gene_features'] + input_list['gene_features_to_compute']

    if config['orthology_table']:
        input_list['orthology_table'] = config['orthology_table']

    if config['ultrametric_species_tree']:
        input_list['ultrametric_species_tree'] = config['ultrametric_species_tree']

    if config['gff']:
        input_list['gff'] = config['gff']

    if config['cafe_results']:
        input_list['cafe_results'] = config['cafe_results']

    if config['focus_species']:
        input_list['focus_species'] = config['focus_species']

    if config['focus_gene_list_files']:
        input_list['focus_gene_list_files'] = config['focus_gene_list_files']

    if config['gene_sets']:
        input_list['gene_sets'] = config['gene_sets']

    if config['orthogroup_sets']:
        input_list['orthogroup_sets'] = config['orthogroup_sets']

    if config['go_universe']:
        input_list['go_universe'] = config['go_universe']

    if config['ontology']:
        input_list['ontology'] = config['ontology']

    return input_list
