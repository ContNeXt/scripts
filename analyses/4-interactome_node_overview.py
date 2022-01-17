# -*- coding: utf-8 -*-

"""Protein-level analysis for ContNeXt: Node overview table of each network"""
# TODO write function docstrings

import os

import matplotlib as mpl
import networkx as nx
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm

from network_utils import create_network_from_edge_file, edge_file_path, load_interactome

mpl.rcParams['figure.dpi'] = 320

# replace here the location of the external data dir if not structured as instructed
data_dir = os.path.join(os.path.expanduser("~"), "contnext_data", "data")

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

# interactome PPI network
interactome = load_interactome(os.path.join(data_dir, "interactome", "interactome_18_01_2021.tsv"))

# load controllability analysis results
with open(os.path.join(data_dir, "controllability_analysis", "interactome_node_classifications.tsv")) as f:
    f.readline()  # disregard header
    node_types = {k: v for k, v in (line.split() for line in f)}

# load housekeeping genes
with open(os.path.join(data_dir, "misc_data", "HK_genes.txt"), 'r') as f:
    hk_genes = [line.split()[0] for line in f.readlines()]


# ### MAKE TABLES ####

def node_qualities(net, betweenness_centrality_dict=None, interactome=False):
    node_qualities_dict = {node: {} for node in net.nodes}
    for node in net.nodes:
        if interactome:
            node_qualities_dict[node]["betweenness centrality"] = round(betweenness_centrality_dict[node], 5)
            node_qualities_dict[node]["betweenness centrality rank"] = None
        node_qualities_dict[node]["degree"] = net.degree(node)
        node_qualities_dict[node]["degree rank"] = None
        if interactome:
            node_qualities_dict[node]["controllability"] = node_types[node]
        node_qualities_dict[node]["housekeeping gene"] = node in hk_genes
    node_qualities_df = pd.DataFrame.from_dict(node_qualities_dict, orient="index")
    if interactome:
        node_qualities_df["betweenness centrality rank"] = node_qualities_df["betweenness centrality"].rank(ascending=False)
    node_qualities_df["degree rank"] = node_qualities_df["degree"].rank(ascending=False)
    return node_qualities_df


# Interactome
os.makedirs(os.path.join(data_dir, "node_properies"), exist_ok=True)
interactome_betweenness_centrality = nx.betweenness_centrality(interactome)
interactome_node_df = node_qualities(interactome, interactome_betweenness_centrality, True)
interactome_node_df.astype(str).to_csv(os.path.join(data_dir, "node_properties", 'interactome_node_properties.tsv'), sep='\t', index=True)
# interactome_node_df = pd.read_csv(os.path.join(data_dir, "node_properties", 'interactome_node_properties.tsv'), sep='\t', index_col=0)


# Normalize betweenness centrality
scaler=MinMaxScaler()
interactome_node_df["betweenness centrality norm"]=scaler.fit_transform(interactome_node_df["betweenness centrality"].to_numpy().reshape(-1, 1))

critical_df = interactome_node_df[interactome_node_df.controllability == "critical"]
ordinary_df = interactome_node_df[interactome_node_df.controllability == "ordinary"]
redundant_df = interactome_node_df[interactome_node_df.controllability == "redundant"]
interactome_overview = pd.DataFrame([[len(critical_df), critical_df["betweenness centrality norm"].mean(),
                                      critical_df["betweenness centrality norm"].median(),
                                      critical_df["betweenness centrality norm"].mode()[0],
                                      critical_df["degree"].mean(), critical_df["degree"].median(),
                                      critical_df["degree"].mode()[0],
                                      "{:.2%}".format(
                                          len(critical_df.loc[critical_df["housekeeping gene"] == True]) / len(
                                              critical_df))],

                                     [len(redundant_df), redundant_df["betweenness centrality norm"].mean(),
                                      redundant_df["betweenness centrality norm"].median(),
                                      redundant_df["betweenness centrality norm"].mode()[0],
                                      redundant_df["degree"].mean(), redundant_df["degree"].median(),
                                      redundant_df["degree"].mode()[0],
                                      "{:.2%}".format(
                                          len(redundant_df.loc[redundant_df["housekeeping gene"] == True]) / len(
                                              redundant_df))],

                                     [len(ordinary_df), ordinary_df["betweenness centrality norm"].mean(),
                                      ordinary_df["betweenness centrality norm"].median(),
                                      ordinary_df["betweenness centrality norm"].mode()[0],
                                      ordinary_df["degree"].mean(), ordinary_df["degree"].median(),
                                      ordinary_df["degree"].mode()[0],
                                      "{:.2%}".format(
                                          len(ordinary_df.loc[ordinary_df["housekeeping gene"] == True]) / len(
                                              ordinary_df))],
                                     ],
                                    columns=["total number", "normalized betweenness centrality mean",
                                             'normalized betweenness centrality median', 'normalized betweenness centrality mode',
                                             "degree mean", "degree median", "degree mode",
                                             "proportion housekeeping gene"],
                                    index=["indispensable", 'dispensable',
                                           'neutral'])

# Tissues, cell types, and cell lines
os.makedirs(os.path.join(data_dir, "node_properies", "tissue"), exist_ok=True)
os.makedirs(os.path.join(data_dir, "node_properies", "cell_type"), exist_ok=True)
os.makedirs(os.path.join(data_dir, "node_properies", "cell_line"), exist_ok=True)

for context, net_dict in zip(['tissue', 'cell_type', 'cell_line'],
                             [tis_network_dict, ct_network_dict, cl_network_dict]):
    for ID, net in tqdm(net_dict.items(), desc=f"iteracting {context} nets"):
        fname = os.path.join(data_dir, "node_properies", context, 'node_properties.tsv')
        node_df = node_qualities(net)
        node_df.astype(str).to_csv(fname, sep='\t', index=True)
