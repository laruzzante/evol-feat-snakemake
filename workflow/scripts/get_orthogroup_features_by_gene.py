import pickle
import statistics

merged_orthogroup_features = snakemake.input.merge_orthogroup_features
orthogroups = pickle.load(open(snakemake.input.orthogroups, 'rb'))
genes = pickle.load(open(snakemake.input.orthogroups, 'rb'))
outfile = snakemake.output[0]
