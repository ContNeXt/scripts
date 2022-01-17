# -*- coding: utf-8 -*-

"""Network-level analysis for ContNeXt: Pairwise co-expression network similarity"""
# TODO write function docstrings

import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from network_utils import create_network_from_edge_file, edge_file_path

mpl.rcParams['figure.dpi'] = 320

# replace here the location of the external data dir if not structured as instructed
data_dir = os.path.join(os.path.expanduser("~"), "contnext_data", "data")

# optional, replace here the desired location of the output figures
figures_dir = os.path.join(data_dir, "figures")
os.makedirs(figures_dir, exist_ok=True)

# ### LOAD DATA ####
# tissue networks
tis_network_dict = {
    ID: create_network_from_edge_file(
        edge_file_path(os.path.join(data_dir, "coexpr_networks", "tissue"), ID), ID)
    for ID in tqdm(os.listdir(os.path.join(data_dir, "coexpr_networks", "tissue")),
                   desc="Creating/loading network objects") if ID != ".DS_Store"
}

# cell type networks
ct_network_dict = {
    ID: create_network_from_edge_file(edge_file_path(os.path.join(data_dir, "coexpr_networks", "cell_type"), ID),
                                      ID)
    for ID in tqdm(os.listdir(os.path.join(data_dir, "coexpr_networks", "cell_type")),
                   desc="Creating/loading network objects") if ID != ".DS_Store"
}

# cell line  networks
cl_network_dict = {
    ID: create_network_from_edge_file(edge_file_path(os.path.join(data_dir, "coexpr_networks", "cell_line"), ID),
                                      ID)
    for ID in tqdm(os.listdir(os.path.join(data_dir, "coexpr_networks", "cell_line")),
                   desc="Creating/loading network objects") if ID != ".DS_Store"
}

# ontology name mappings
with open(os.path.join(data_dir, "mappings", "uberon_name_mappings.json"), 'r') as f:
    uberon_name_mappings = json.load(f)
uberon_id_mappings = {v: k for k, v in uberon_name_mappings.items()}
with open(os.path.join(data_dir, "mappings", "CL_name_mappings.json"), 'r') as f:
    CL_name_mappings = json.load(f)
CL_id_mappings = {v: k for k, v in CL_name_mappings.items()}
with open(os.path.join(data_dir, "mappings", "CLO_name_mappings.json"), 'r') as f:
    CLO_name_mappings = json.load(f)
CLO_id_mappings = {v: k for k, v in CLO_name_mappings.items()}

# ### PAIRWISE NET EDGE COMPARISONS ####
# Tissues
tissues = [uberon_name_mappings[ID] for ID in tis_network_dict]
similarity_tissues = pd.DataFrame(index=tissues, columns=tissues)

for i, row in tqdm(similarity_tissues.iterrows(), total=len(similarity_tissues), desc="calculating similarity by row"):
    net1 = tis_network_dict[uberon_id_mappings[i]]
    for tis in similarity_tissues.columns:
        if pd.notnull(similarity_tissues.loc[tis, i]):
            continue
        net2 = tis_network_dict[uberon_id_mappings[tis]]
        if net1 == net2:
            intersection = len(net1.edges())
        else:
            intersection_net = nx.intersection(net1.subgraph(net2.nodes), net2.subgraph(net1.nodes))
            intersection = len(intersection_net.edges())
        row[tis] = intersection
similarity_tissues.to_csv(os.path.join(data_dir, "misc_data", 'pairwise_similarity_tissues.tsv'), sep='\t', index=True)
# if already done, load it
# similarity_tissues = pd.read_csv(os.path.join(data_dir, "misc_data", 'pairwise_similarity_tissues.tsv'), sep='\t', index_col=0)

# Cell types
cell_types = [CL_name_mappings[ID] for ID in ct_network_dict]
similarity_cell_types = pd.DataFrame(index=cell_types, columns=cell_types)

for i, row in tqdm(similarity_cell_types.iterrows(), total=len(similarity_cell_types),
                   desc="calculating similarity by row"):
    net1 = ct_network_dict[CL_id_mappings[i]]
    for ct in similarity_cell_types.columns:
        if pd.notnull(similarity_cell_types.loc[ct, i]):
            continue
        net2 = ct_network_dict[CL_id_mappings[ct]]
        if net1 == net2:
            intersection = len(net1.edges())
        else:
            intersection_net = nx.intersection(net1.subgraph(net2.nodes), net2.subgraph(net1.nodes))
            intersection = len(intersection_net.edges())
        row[ct] = intersection
similarity_cell_types = similarity_cell_types.apply(pd.to_numeric)
similarity_cell_types.to_csv(os.path.join(data_dir, "misc_data", 'pairwise_similarity_cell_types.tsv'), sep='\t', index=True)

# if alread done, load it
# similarity_cell_types = pd.read_csv(os.path.join(data_dir, "misc_data", 'pairwise_similarity_cell_types.tsv'), sep='\t', index_col=0)


# Cell lines
cell_lines = [CLO_name_mappings[ID] for ID in cl_network_dict]
similarity_cell_lines = pd.DataFrame(index=cell_lines, columns=cell_lines)

for i, row in tqdm(similarity_cell_lines.iterrows(), total=len(similarity_cell_lines),
                   desc="calculating similarity by row"):
    net1 = cl_network_dict[CLO_id_mappings[i]]
    for cl in similarity_cell_lines.columns:
        if pd.notnull(similarity_cell_lines.loc[cl, i]):
            continue
        net2 = cl_network_dict[CLO_id_mappings[cl]]
        if net1 == net2:
            intersection = len(net1.edges())
        else:
            intersection_net = nx.intersection(net1.subgraph(net2.nodes), net2.subgraph(net1.nodes))
            intersection = len(intersection_net.edges())
        row[cl] = intersection
similarity_cell_lines = similarity_cell_lines.apply(pd.to_numeric)
similarity_cell_lines.to_csv(os.path.join(data_dir, "misc_data", 'pairwise_similarity_cell_lines.tsv'), sep='\t', index=True)
# if alread done, load it
# similarity_cell_lines = pd.read_csv(os.path.join(data_dir, "misc_data", 'pairwise_similarity_cell_lines.tsv'), sep='\t', index_col=0)


# ### FIGURES ####
# Similarity matrices as heatmaps

blue_color_palette = sns.color_palette("Blues", as_cmap=True)


def plot_df(df, context, figname=None):
    mpl.rcParams['xtick.labelsize'] = 'medium'
    mpl.rcParams['ytick.labelsize'] = 'medium'
    g = sns.clustermap(df, row_cluster=False, col_cluster=False, cmap=blue_color_palette, figsize=(20, 15),
                       cbar_pos=None,
                       yticklabels=df.index,
                       xticklabels=df.index,
                       )
    g.ax_col_dendrogram.set_visible(False)
    g.ax_row_dendrogram.set_visible(False)
    plt.title(f"Pairwise similarity matrix for {context} co-expression networks", fontsize=30, weight="bold")
    if figname:
        plt.savefig(figname, bbox_inches='tight', dpi=320)
    return g


# remove diagonal
for i in range(len(similarity_tissues.index)):
    similarity_tissues.iloc[i, i] = 0

for i in range(len(similarity_cell_types.index)):
    similarity_cell_types.iloc[i, i] = 0

for i in range(len(similarity_cell_lines.index)):
    similarity_cell_lines.iloc[i, i] = 0

plot_df(similarity_tissues, "tissue", os.path.join(figures_dir, "pairwise_tissue_similarity.png"))
plot_df(similarity_cell_types, "cell type", os.path.join(figures_dir, "pairwise_celltype_similarity.png"))
plot_df(similarity_cell_lines, "cell line", os.path.join(figures_dir, "pairwise_cellline_similarity.png"))
