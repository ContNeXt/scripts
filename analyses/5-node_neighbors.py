# -*- coding: utf-8 -*-

"""Protein-level analysis for ContNeXt: node neighbors"""

import json
import os
from collections import defaultdict, Counter

from tqdm import tqdm

from network_utils import create_network_from_edge_file, edge_file_path

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


# ### CREATE NODE NEIGHBOR TABLE ####
for context, net_dict in zip(['tissue', 'cell-type', 'cell-line'], [tis_network_dict, ct_network_dict, cl_network_dict]):
    node_dict = defaultdict(list)
    for net in net_dict.values(): 
        for node in net.nodes():
            for neighbor in net.neighbors(node):
                node_dict.append(neighbor)
    # This is just reducing the list to a dict where you have the # of times that a node appears
    node_dict = {
        node: Counter(neighbors).most_common()
        for node, neighbors in node_dict.items()
    }
    # You can apply a threshold here if you want and take the top 5 neighbors or something
    with open(os.path.join(data_dir, "misc_data", f'{context}_neighbors.json'), 'w') as f:
        json.dump(node_dict, f)
