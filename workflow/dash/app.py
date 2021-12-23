import dash
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import numpy as np
import somoclu
import time
import subprocess
import plotly.graph_objs as go
# import statistics
from collections import defaultdict
from sklearn import preprocessing
from sklearn.cluster import KMeans
# from sklearn.cluster import DBSCAN  # could be OPTICS, SpectralClustering, ..., look at sklearn manual for available options

PORT = 8000
ADDRESS = '127.0.0.1'

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
data = pd.read_csv('../output/merged_orthogroup_features.tsv', sep='\t')
gene_data = pd.read_csv('../output/merged_orthogroup_features.tsv', sep='\t')

# Processing data, removing rows containing NAs, necessary to continue computation as somoclu requires a numpy clean array.
data = data.dropna()

# Reindexing without NA rows to be able to associate som bmus rows to original data rows.
data = data.reset_index(drop=True)

# Removing first column that contains OG ids.
metrics = data.iloc[:, 1:]

# Generating metrics names list for dropdown menu
metrics_names = list(metrics.columns.values)
metrics_names_dict = []
z = 0
for metric in metrics_names:
    metrics_names_dict.append({'label': metric, 'value': z})
    z += 1

# Generating metrics names list for dropdown menu
#gene_ids = list(gene_data['gene_id'])
#gene_ids_dict = []
#for gene_id in gene_ids:
#    gene_ids_dict.append({'label': gene_id, 'value': gene_id})

# Scaling metrics
scaled_metrics = preprocessing.scale(metrics, axis=0)  # axis = 0 to standardize each feature (column) independently.

"""
Code to verify that the scale function does indeed put mean to 0 and stdev to 1 per column.
for i in range(np.size(scaled_metrics, 1)):
    print(i)
    print(statistics.mean(scaled_metrics[0:, i]))
    print(statistics.stdev(scaled_metrics[0:, i]))
"""

# Converting to numpy array, required for somoclu
scaled_metrics = np.float32(scaled_metrics)

# Setting SOM grid size
n_rows, n_cols = 25, 25

# Computing parallel computation of SOM
som = somoclu.Somoclu(n_cols, n_rows, maptype="toroid", compactsupport=False, initialization=None)
som.train(scaled_metrics)

# Computing superclusters of SOM cells
start_time = time.time()
som.cluster()  # by default uses KMeans algorithm with 8 clusters.
# algorithm = DBSCAN()  # example of using a different clustering algoritm
# som.cluster(algorithm=algorithm)
print("--- %s seconds ---" % (time.time() - start_time))

# Computing map from som bmus coordinates to data row indices, in order to find OGs in raw data by clicking heatmap tiles (i.e. getting som bmus coordinates).
row_ind_2_coords = {}
coords_2_row_ind = defaultdict(list)
for i in data.index.values:
    row_ind_2_coords[i] = som.bmus[i]
    coords_2_row_ind[som.bmus[i][0], som.bmus[i][1]].append(i)

# Computing dicionary of codebook vectors matrices, for each Metric separately, from som codebook lists. With this, we can plot the heatmap using the correct required matrix structure.
codebook_matrices = defaultdict(list)
for z in range(len(metrics.columns)):
    matrix = []
    for i in range(n_rows):
        row = []
        for j in range(n_cols):
            row.append(som.codebook[i][j][z])
        matrix.append(row)
    codebook_matrices[z] = matrix
width, height = 500, 500


def generate_table(dataframe, max_rows=5):
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +
        # Body
        [html.Tr([
            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
        ]) for i in range(min(len(dataframe), max_rows))]
    )


def generate_boxplots(dataframe):
    boxplots = []
    for metric in dataframe.columns:
        if metric != 'orthogroup':
            boxplots.append(go.Box(
                x=preprocessing.scale(dataframe[metric]),
                name=metric,
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
    ], className="six columns", style={'diplay': 'inline-block'}),
    html.Div([
        dcc.Graph(
            id='feature-heatmap',
        ),
        html.Label('Select feature:'),
        dcc.Dropdown(
            id='feature-dropdown',
            options=metrics_names_dict,
            value=0
        ),
        html.H2(' '),
        html.Label('Select gene:'),
        dcc.Dropdown(
            id='gene-dropdown',
            #options=gene_ids_dict,
            value=None
        ),
        generate_table(data),
        dash_table.DataTable(id='go-table', columns=[])
    ], className="six columns")
], className="row")


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


@app.callback(
    output=dash.dependencies.Output('feature-boxplots', 'figure'),
    inputs=[dash.dependencies.Input('cluster-heatmap', 'clickData')])
def update_boxplots(clickData):
        print(clickData)
        row_indices = coords_2_row_ind[clickData['points'][0]['y'], clickData['points'][0]['x']]
        selected_data = data.iloc[row_indices]
        print(data.iloc[row_indices])
        #selected_gene_data = gene_data.loc[gene_data['pub_og_id'].isin(selected_data['pub_og_id'])]
        print(selected_gene_data)
        #genesOfInterest = selected_gene_data['gene_id']
        #f = open('genesOfInterest.tsv', 'w')
        f.write('\n'.join(genesOfInterest))
        f.close()
        #subprocess.call(['Rscript', '~/evol-feat/go_enrichment_scripts/topGO_geneList.R'])
        generate_table(selected_gene_data)
        return{
            'data': generate_boxplots(selected_data)
        }


@app.callback(
    dash.dependencies.Output('feature-heatmap', 'figure'),
    [dash.dependencies.Input('feature-dropdown', 'value')])
def update_feature_heatmap(metric):
    return {
        'data': [{
            'z': codebook_matrices[metric],
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
