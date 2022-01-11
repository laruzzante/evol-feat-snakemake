import dash
import dash_table
# import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import pandas as pd
import numpy as np
import somoclu
import time
import subprocess
import pickle
import plotly.graph_objs as go
# import statistics
from collections import defaultdict
from sklearn import preprocessing
from sklearn.cluster import KMeans
# from sklearn.cluster import DBSCAN  # could be OPTICS, SpectralClustering, ..., look at sklearn manual for available options

PORT = 8000
ADDRESS = '127.0.0.1'

PAGE_SIZE = 10


external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
orthogroups_data = pd.read_csv('../output/merged_orthogroup_features.tsv', sep='\t')
orthogroups_dict = pickle.load(open('../output/.orthogroups.pickle', 'rb'))
genes_dict = pickle.load(open('../output/.genes.pickle', 'rb'))
species_dict = pickle.load(open('../output/.species.pickle', 'rb'))


# Processing data, removing rows containing NAs, necessary to continue computation as somoclu requires a numpy clean array.
orthogroups_data = orthogroups_data.dropna()

# Reindexing without NA rows to be able to associate som bmus rows to original data rows.
orthogroups_data = orthogroups_data.reset_index(drop=True)

# Removing first column that contains OG ids.
features = orthogroups_data.iloc[:, 1:]

# Generating features names list for dropdown menu
feature_names = list(features.columns.values)
feature_names_list = []
z = 0
for feature in sorted(feature_names):
    feature_names_list.append({'label': feature, 'value': z})
    z += 1

# # Generating gene ids names list for dropdown menu
# gene_ids = list(genes_dict.keys())
# gene_ids_list = []
# for gene_id in genes_dict.keys():
#    gene_ids_list.append({'label': gene_id, 'value': gene_id})

# Generating species ids names list for dropdown menu
species_ids = list(species_dict.keys())
species_ids_list = []
for species_id in sorted(species_dict.keys()):
   species_ids_list.append({'label': species_id, 'value': species_id})

# Scaling features
scaled_features = preprocessing.scale(features, axis=0)  # axis = 0 to standardize each feature (column) independently.

"""
Code to verify that the scale function does indeed put mean to 0 and stdev to 1 per column.
for i in range(np.size(scaled_features, 1)):
    print(i)
    print(statistics.mean(scaled_features[0:, i]))
    print(statistics.stdev(scaled_features[0:, i]))
"""

# Converting to numpy array, required for somoclu
scaled_features = np.float32(scaled_features)

# Setting SOM grid size
n_rows, n_cols = 25, 25

# Computing parallel computation of SOM
som = somoclu.Somoclu(n_cols, n_rows, maptype="toroid", compactsupport=False, initialization=None)
som.train(scaled_features)

# Computing superclusters of SOM cells
start_time = time.time()
som.cluster()  # by default uses KMeans algorithm with 8 clusters.
# algorithm = DBSCAN()  # example of using a different clustering algoritm
# som.cluster(algorithm=algorithm)
print("--- %s seconds ---" % (time.time() - start_time))

# Computing map from som bmus coordinates to data row indices, in order to find OGs in raw data by clicking heatmap tiles (i.e. getting som bmus coordinates).
row_ind_2_coords = {}
coords_2_row_ind = defaultdict(list)
for i in orthogroups_data.index.values:
    row_ind_2_coords[i] = som.bmus[i]
    coords_2_row_ind[som.bmus[i][0], som.bmus[i][1]].append(i)

# Computing dicionary of codebook vectors matrices, for each Feature separately, from som codebook lists. With this, we can plot the heatmap using the correct required matrix structure.
codebook_matrices = defaultdict(list)
for z in range(len(features.columns)):
    matrix = []
    for i in range(n_rows):
        row = []
        for j in range(n_cols):
            row.append(som.codebook[i][j][z])
        matrix.append(row)
    codebook_matrices[z] = matrix
width, height = 500, 500


# def generate_table(dataframe, max_rows=5):
#     return html.Table(
#         # Header
#         [html.Tr([html.Th(col) for col in dataframe.columns])] +
#         # Body
#         [html.Tr([
#             html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
#         ]) for i in range(min(len(dataframe), max_rows))]
#     )


def generate_boxplots(dataframe):
    boxplots = []
    for feature in dataframe.columns:
        if feature != 'orthogroup':
            boxplots.append(go.Box(
                x=preprocessing.scale(dataframe[feature]),
                name=feature,
                jitter=0.3,
                pointpos=0,
                boxpoints='all',
                marker=dict(
                    color='rgb(7,40,89)'),
                line=dict(
                    color='rgb(209,44,209)')
            ))
    return boxplots


app.layout = html.Div([
    html.Div([
        dcc.Graph(
            id='cluster-heatmap',
        ),
        html.Label('KMeans N:'),
        dcc.Slider(
            id='clusters-slider',
            min=2,
            max=15,
            step=1,
            value=8,
            marks={i: '{}'.format(i) for i in range(16)},
            included=False,
        ),
        html.H2(" "),
        dcc.Graph(
            id='feature-boxplots'
        )
    ], className="six columns", style={'display': 'inline-block'}),

    html.Div([
        dcc.Graph(
            id='feature-heatmap',
        ),

        html.Label('Select feature:'),
        dcc.Dropdown(
            id='feature-dropdown',
            options=feature_names_list
        ),

        html.H2(' '),

        html.Label('Select species:'),
        dcc.Dropdown(
            id='species-dropdown',
            options=species_ids_list
        ),


        # html.Label('Select gene:'),
        # dcc.Dropdown(
        #     id='gene-dropdown',
        #     options=gene_ids_list,
        # ),

        # generate_table(round(orthogroups_data, 2)),

        dash_table.DataTable(
            id='table-sorting-filtering',
            columns=[
                {'name': i, 'id': i, 'deletable': True} for i in orthogroups_data.columns
            ],
            page_current= 0,
            page_size=PAGE_SIZE,
            page_action='custom',

            filter_action='custom',
            filter_query='',

            sort_action='custom',
            sort_mode='multi',
            sort_by=[]
        ),


        dash_table.DataTable(id='go-table', columns=[])
    ], className="six columns")
], className="row")



operators = [['ge ', '>='],
             ['le ', '<='],
             ['lt ', '<'],
             ['gt ', '>'],
             ['ne ', '!='],
             ['eq ', '='],
             ['contains ']]

def split_filter_part(filter_part):
    for operator_type in operators:
        for operator in operator_type:
            if operator in filter_part:
                name_part, value_part = filter_part.split(operator, 1)
                name = name_part[name_part.find('{') + 1: name_part.rfind('}')]

                value_part = value_part.strip()
                v0 = value_part[0]
                if (v0 == value_part[-1] and v0 in ("'", '"', '`')):
                    value = value_part[1: -1].replace('\\' + v0, v0)
                else:
                    try:
                        value = float(value_part)
                    except ValueError:
                        value = value_part

                # word operators need spaces after them in the filter string,
                # but we don't want these later
                return name, operator_type[0].strip(), value

    return [None] * 3



@app.callback(
    Output('table-sorting-filtering', 'data'),
    Input('table-sorting-filtering', "page_current"),
    Input('table-sorting-filtering', "page_size"),
    Input('table-sorting-filtering', 'sort_by'),
    Input('table-sorting-filtering', 'filter_query'))
def update_table(page_current, page_size, sort_by, filter):
    filtering_expressions = filter.split(' && ')
    dff = round(orthogroups_data, 2)
    for filter_part in filtering_expressions:
        col_name, operator, filter_value = split_filter_part(filter_part)

        if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
            # these operators match pandas series operator method names
            dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
        elif operator == 'contains':
            dff = dff.loc[dff[col_name].str.contains(filter_value)]

    if len(sort_by):
        dff = dff.sort_values(
            [col['column_id'] for col in sort_by],
            ascending=[
                col['direction'] == 'asc'
                for col in sort_by
            ],
            inplace=False
        )

    page = page_current
    size = page_size
    return dff.iloc[page * size: (page + 1) * size].to_dict('records')



@app.callback(
    dash.dependencies.Output('cluster-heatmap', 'figure'),
    [dash.dependencies.Input('clusters-slider', 'value')])
def update_cluster_heatmap(n_clusters):
    som.cluster(algorithm=KMeans(n_clusters))
    return {
        'data': [{
            'z': som.clusters,
            'type': 'heatmap',
        }],
        'layout': {
            'title': 'Evol-Feat SOM by Cluster',
            'width': width,
            'height': height
        }
    }


# @app.callback(
#     output=dash.dependencies.Output('feature-boxplots', 'figure'),
#     inputs=[dash.dependencies.Input('cluster-heatmap', 'clickData')])
# def update_boxplots(clickData):
#         print(clickData)
#         row_indices = coords_2_row_ind[clickData['points'][0]['y'], clickData['points'][0]['x']]
#         selected_data = orthogroups_data.iloc[row_indices]
#         print(orthogroups_data.iloc[row_indices])
#         #selected_gene_data = gene_data.loc[gene_data['pub_og_id'].isin(selected_data['pub_og_id'])]
#         #print(selected_gene_data)
#         #genesOfInterest = selected_gene_data['gene_id']
#         #f = open('genesOfInterest.tsv', 'w')
#         f.write('\n'.join(genesOfInterest))
#         f.close()
#         #subprocess.call(['Rscript', '~/evol-feat/go_enrichment_scripts/topGO_geneList.R'])
#         generate_table(selected_gene_data)
#         return{
#             'data': generate_boxplots(selected_data)
#         }


@app.callback(
    dash.dependencies.Output('feature-heatmap', 'figure'),
    [dash.dependencies.Input('feature-dropdown', 'value')])
def update_feature_heatmap(feature):
    return {
        'data': [{
            'z': codebook_matrices[feature],
            'type': 'heatmap'
        }],
        'layout': {
            'title': 'Evol-Feat SOM by Feature',
            'width': width,
            'height': height
        },
    }

# @app.callback(
#     dash.dependencies.Output(‘go-table’, ‘data’),
#     [dash.dependencies.Input(‘gene-dropdown’, ‘value’)])
# def update_go_table(gene_id):
#     return {
#         columns=[{“name”: i, “id”: i} for i in gene_data.columns],
#         data=df.to_dict(‘records’)
#     }


if __name__ == '__main__':
    app.run_server(
        port=PORT,
        host=ADDRESS)
