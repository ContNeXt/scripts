# -*- coding: utf-8 -*-

"""Network-level analysis for ContNeXt:
Calculated the most common edges for co-expression networks"""

import json
import os

from tqdm import tqdm

from network_utils import create_network_from_edge_file, edge_file_path, count_in_freq, most_common

# replace here the location of the external data dir if not structured as instructed
data_dir = os.path.join(os.path.expanduser("~"), "contnext_data", "data")


# ### TISSUES ####
print("---Tissues---")
tis_network_dict = {
    ID: create_network_from_edge_file(
        edge_file_path(os.path.join(data_dir, "coexpr_networks", "tissue"), ID), ID)
    for ID in tqdm(os.listdir(os.path.join(data_dir, "coexpr_networks", "tissue")),
                   desc="Creating/loading network objects") if ID != ".DS_Store"
}

print("calculating most common edges...")
tis_most_common_edges = most_common(tis_network_dict.values(), comparison="edges")


print("...saving 100,000 most common edges to json...")

with open(os.path.join(data_dir, "misc_data", "tis_100000most_common_edges.json"), 'w') as f:
    json.dump(tis_most_common_edges[:100000], f)

print("...counting unique edges...")
count = 0
for e, freq in tis_most_common_edges:
    if count_in_freq(freq) == 1:
        count += 1
print(f"\tthere are {count} unique edges")

print("...complete")


# ### CELL TYPES ####
print("---Cell types---")
ct_network_dict = {
    ID: create_network_from_edge_file(edge_file_path(os.path.join(data_dir, "coexpr_networks", "cell_type"), ID),
                                      ID)
    for ID in tqdm(os.listdir(os.path.join(data_dir, "coexpr_networks", "cell_type")),
                   desc="Creating/loading network objects") if ID != ".DS_Store"
}

print("calculating most common edges...")
ct_most_common_edges = most_common(ct_network_dict.values(), comparison="edges")


print("...saving 100,000 most common edges to json...")

with open(os.path.join(data_dir, "misc_data", "ct_100000most_common_edges.json"), 'w') as f:
    json.dump(ct_most_common_edges[:100000], f)

print("...counting unique edges...")
count = 0
for e, freq in ct_most_common_edges:
    if count_in_freq(freq) == 1:
        count += 1
print(f"\tthere are {count} unique edges")

print("...complete")


# ### CELL TYPES ####
print("---Cell types---")
cl_network_dict = {
    ID: create_network_from_edge_file(edge_file_path(os.path.join(data_dir, "coexpr_networks", "cell_line"), ID),
                                      ID)
    for ID in tqdm(os.listdir(os.path.join(data_dir, "coexpr_networks", "cell_line")),
                   desc="Creating/loading network objects") if ID != ".DS_Store"
}

print("calculating most common edges...")
cl_most_common_edges = most_common(cl_network_dict.values(), comparison="edges")


print("...saving 100,000 most common edges to json...")

with open(os.path.join(data_dir, "misc_data", "cl_100000most_common_edges.json"), 'w') as f:
    json.dump(cl_most_common_edges[:100000], f)

print("...counting unique edges...")
count = 0
for e, freq in cl_most_common_edges:
    if count_in_freq(freq) == 1:
        count += 1
print(f"\tthere are {count} unique edges")

print("...complete")