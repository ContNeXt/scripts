# -*- coding: utf-8 -*-

"""Network-level analysis for ContNeXt:
Co-expression patterns correspond with protein-protein interactions 'more than random' """

import json
import os
import random
import sys

import networkx as nx
import pandas as pd
from tqdm import tqdm
from xswap.permute import permute_edge_list
from xswap.preprocessing import map_str_edges

from network_utils import create_network_from_edge_file, edge_file_path, load_interactome

num_permutations = int(sys.argv[1])

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

# name mappings
with open(os.path.join(data_dir, "mappings", "uberon_name_mappings.json"), 'r') as f:
    uberon_name_mappings = json.load(f)
with open(os.path.join(data_dir, "mappings", "CL_name_mappings.json"), 'r') as f:
    CL_name_mappings = json.load(f)
with open(os.path.join(data_dir, "mappings", "CLO_name_mappings.json"), 'r') as f:
    CLO_name_mappings = json.load(f)
with open(os.path.join(data_dir, "mappings", "hgnc_name_mappings.json"), 'r') as f:
    hgnc_mappings = json.load(f)


# ### ANALYSIS ####
# permute network
def permute_network(network) -> nx.Graph:
    """Permute network.

    :param network: nx network object
    :return: permuted network
    """
    # Read edge list
    edge_list = network.edges()
    # Get mapping since the edge list contains now integers
    edge_list_integers, node_mapping, _ = map_str_edges(edge_list, bipartite=False)
    # Permute
    permuted_edges, stats = permute_edge_list(edge_list_integers, seed=random.randint(1, 5000))
    # Reverse mapping dictionary
    node_mapping = {
        v: k
        for k, v in node_mapping.items()
    }
    new_edges = [
        (node_mapping[source], node_mapping[target])
        for index, (source, target) in enumerate(permuted_edges)
    ]
    new_net = nx.Graph()
    new_net.add_edges_from(new_edges)
    return new_net


# create table for edge overlap of each permuted net to the interactome
f = open(os.path.join(data_dir, "misc_data", "permuted_nets_similarity_to_interactome.tsv"), "w")
for context, net_dict in zip(['tissue', 'cell type', 'cell line'],
                             [tis_network_dict, ct_network_dict, cl_network_dict]):
    for ID, net in tqdm(net_dict.items(), desc=f"iterating {context} nets"):
        f.write(f"{context}\t{ID}\t")
        for x in range(num_permutations):
            count = 0
            permuted_net = permute_network(net)
            intersection_net = nx.intersection(permuted_net.subgraph(interactome.nodes),
                                               interactome.subgraph(permuted_net.nodes))
            intersection = len(intersection_net.edges())
            f.write(f"{intersection}\t")
        f.write("\n")
f.close()

# load results
results = pd.read_csv(os.path.join(data_dir, "misc_data", "permuted_nets_similarity_to_interactome.tsv"), header=None, sep="\t",
                      dtype={1: str})
results = results.drop(102, axis='columns')
results.columns = ["Context", "ID", "Overlap after permutations (1...)", *list(range(2, 101))]

# add unpermuted true overlap
true_overlap_with_interactome = {}
for ID, net in tqdm(tis_network_dict.items()):
    intersection_net = nx.intersection(net.subgraph(interactome.nodes), interactome.subgraph(net.nodes))
    intersection = len(intersection_net.edges())
    true_overlap_with_interactome[ID] = intersection
for ID, net in tqdm(ct_network_dict.items()):
    intersection_net = nx.intersection(net.subgraph(interactome.nodes), interactome.subgraph(net.nodes))
    intersection = len(intersection_net.edges())
    true_overlap_with_interactome[ID] = intersection
for ID, net in tqdm(cl_network_dict.items()):
    intersection_net = nx.intersection(net.subgraph(interactome.nodes), interactome.subgraph(net.nodes))
    intersection = len(intersection_net.edges())
    true_overlap_with_interactome[ID] = intersection
new_col = [true_overlap_with_interactome[ID] for ID in results["ID"]]
results.insert(loc=2, column='Overlap', value=new_col)

# add extra info
avg_overlap = results.iloc[:, 3:].mean(axis=1)
median_overlap = results.iloc[:, 3:].median(axis=1)
results.insert(loc=3, column="Mean overlap after permutations", value=avg_overlap)
results.insert(loc=4, column="Median overlap after permutations", value=median_overlap)

names = []
for ID in results.loc[results["Context"] == "tissue"]["ID"]:
    names.append(uberon_name_mappings[str(ID)])
for ID in results.loc[results["Context"] == "cell type"]["ID"]:
    names.append(CL_name_mappings[str(ID)])
for ID in results.loc[results["Context"] == "cell line"]["ID"]:
    names.append(CLO_name_mappings[str(ID)])
results.insert(loc=2, column="Name", value=names)

net_size = {
    net_ID: len(net.nodes()) for net_ID, net in tis_network_dict.items()
}
for net_ID, net in cl_network_dict.items():
    net_size[net_ID] = len(net.nodes())
sizes = [net_size[ID] for ID in results["ID"]]
results.insert(loc=3, column="Network size", value=sizes)

factor_diff = []
for i, row in results.iterrows():
    true_overlap = row["Overlap"]
    avg_perm_overlap = row["Mean overlap after permutations"]
    factor_diff.append(true_overlap / avg_perm_overlap)
results.insert(loc=4, column="difference", value=factor_diff)

# save final results
results.to_csv(os.path.join(data_dir, "misc_data", 'interactome-net-similarity_permuted.tsv'), sep='\t', index=False)
