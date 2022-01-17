# -*- coding: utf-8 -*-

"""Main protein-level analysis for ContNeXt"""

import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from tqdm import tqdm

from network_utils import create_network_from_edge_file, edge_file_path, load_interactome, most_common, least_common, \
    identify_hub_genes, venn_d, above_cutoff, below_cutoff, get_genes_in_overlap, get_top_x_genes_in_overlap, \
    which_nets_have

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

# interactome PPI network
interactome = load_interactome(os.path.join(data_dir, "interactome", "interactome_18_01_2021.tsv"))

# ontology name mappings
with open(os.path.join(data_dir, "mappings", "uberon_name_mappings.json"), 'r') as f:
    uberon_name_mappings = json.load(f)
with open(os.path.join(data_dir, "mappings", "CL_name_mappings.json"), 'r') as f:
    CL_name_mappings = json.load(f)
with open(os.path.join(data_dir, "mappings", "CLO_name_mappings.json"), 'r') as f:
    CLO_name_mappings = json.load(f)

# load indispensable interactome nodes
with open(os.path.join(data_dir, "controllability_analysis", "interactome_indispensable_nodes.txt"), 'r') as f:
    indispensable_nodes = [line.strip() for line in f.readlines()]

# load housekeeping genes
with open(os.path.join(data_dir, "misc_data", "HK_genes.txt"), 'r') as f:
    hk_genes = [line.split()[0] for line in f.readlines()]

# ### PREP DATA FOR ANALYSIS ####

most_common_nodes_in_tis = most_common(tis_network_dict.values(), comparison="nodes")
least_common_nodes_in_tis = least_common(tis_network_dict.values(), comparison="nodes")

most_common_nodes_in_ct = most_common(ct_network_dict.values(), comparison="nodes")
least_common_nodes_in_ct = least_common(ct_network_dict.values(), comparison="nodes")

most_common_nodes_in_cl = most_common(cl_network_dict.values(), comparison="nodes")
least_common_nodes_in_cl = least_common(cl_network_dict.values(), comparison="nodes")

all_tissue_nodes = above_cutoff(most_common_nodes_in_tis, 0)
all_celltype_nodes = above_cutoff(most_common_nodes_in_ct, 0)
all_cellline_nodes = above_cutoff(most_common_nodes_in_cl, 0)

all_nodes = set(all_tissue_nodes).union(set(all_celltype_nodes)).union(set(all_cellline_nodes))

# Identify the most connected nodes in interactome
pathway_hubgenes = identify_hub_genes(interactome)


# ### MOST PREVALANT PROTEINS ####
# Tissues
in_all_tis = above_cutoff(most_common_nodes_in_tis, 46)
print(f"there are {len(in_all_tis)} proteins in all tissue networks")
unique_tis = below_cutoff(least_common_nodes_in_tis, 1)
print(f"there are {len(unique_tis)} protein unique to one tissue network")
print("indispensable interactome nodes that are in all tissue networks:",
      [node for node in in_all_tis if node in indispensable_nodes])
print(
    f"There are {len([node for node in in_all_tis if node in hk_genes])} housekeeping genes that are in all tissue networks:",
    [node for node in in_all_tis if node in hk_genes])

# Cell types
in_all_ct = above_cutoff(most_common_nodes_in_ct, 30)
print(f"there are {len(in_all_ct)} proteins in all cell type networks")
unique_ct = below_cutoff(least_common_nodes_in_ct, 1)
print(f"there is {len(unique_ct)} protein unique to one cell type network:",  # only 1
      unique_ct)
print(
    f"the {CL_name_mappings[which_nets_have(unique_ct[0], 'node', ct_network_dict)[0]]} network has the {unique_ct[0]} protein")
print(f"{unique_ct[0]} is {'' if unique_ct[0] in indispensable_nodes else 'not'} indispensable")
print("indispensable interactome nodes that are in all cell type networks:",
      ", ".join(sorted([node for node in in_all_ct if node in indispensable_nodes])))
print(
    f"There are {len([node for node in in_all_ct if node in hk_genes])} housekeeping genes that are in all cell type networks")

# Cell lines
in_all_cl = above_cutoff(most_common_nodes_in_cl, 22)
print(f"there are {len(in_all_cl)} proteins in all cell line networks")
unique_cl = below_cutoff(least_common_nodes_in_cl, 1)
print(f"there is {len(unique_cl)} protein unique to one cell type network")
print(
    f"{len([node for node in unique_cl if node in indispensable_nodes])} of the unique cell line nodes are indispensable")
print("indispensable interactome nodes that are in all cell line networks:",
      ", ".join(sorted([node for node in in_all_cl if node in indispensable_nodes])))
print(
    f"There are {len([node for node in in_all_cl if node in hk_genes])} housekeeping genes that are in all cell line networks")

# ### FIGURES ####
# Comparing Interactome vs coexp nets
# Tissues
interactome_vs_tis_overlap = set([n for n, _ in pathway_hubgenes]).intersection(
    set([n for n, _ in most_common_nodes_in_tis]))
print(
    f"The interactome and all the combination of all tissue co-expression network have {len(interactome_vs_tis_overlap)} nodes in common"
    f"\n({len(interactome_vs_tis_overlap) / len(pathway_hubgenes):.2%} of all nodes in the interactome)"
    f"\n({len(interactome_vs_tis_overlap) / len(most_common_nodes_in_tis):.2%} of all nodes in tissue co-expression networks)")
# Cell types
interactome_vs_ct_overlap = set([n for n, _ in pathway_hubgenes]).intersection(
    set([n for n, _ in most_common_nodes_in_ct]))
print(
    f"The interactome and all the combination of all cell type co-expression network have {len(interactome_vs_ct_overlap)} nodes in common"
    f"\n({len(interactome_vs_ct_overlap) / len(pathway_hubgenes):.2%} of all nodes in the interactome)"
    f"\n({len(interactome_vs_ct_overlap) / len(most_common_nodes_in_ct):.2%} of all nodes in cell type co-expression networks)")
# Cell lines
interactome_vs_cl_overlap = set([n for n, _ in pathway_hubgenes]).intersection(
    set([n for n, _ in most_common_nodes_in_cl]))
print(
    f"The interactome and all the combination of all cell line co-expression network have {len(interactome_vs_cl_overlap)} nodes in common"
    f"\n({len(interactome_vs_cl_overlap) / len(pathway_hubgenes):.2%} of all nodes in the interactome)"
    f"\n({len(interactome_vs_cl_overlap) / len(most_common_nodes_in_cl):.2%} of all nodes in cell line co-expression networks)")

# print(interactome_vs_tis_overlap == interactome_vs_ct_overlap)
# print(interactome_vs_tis_overlap == interactome_vs_cl_overlap)
# print(interactome_vs_tis_overlap-interactome_vs_cl_overlap)

venn_d(set([n for n, _ in pathway_hubgenes]), set(all_nodes),
       "All interactome proteins", "All co-expression network proteins", "All Interactome vs all context proteins",
       output_filename=os.path.join(figures_dir, "all_interactome_vs_co-exp-net_proteins.png"))

# Comparing most common proteins accross contexts to most connected proteins in interactome
# Tissues
cutoffs1 = [46, 45, 44, 43]
top_x = [10, 50, 150, 400]
tis_prom_genes = []
tis_most_prom_genes_in_int_overlap = []
tis_top_pathway_hubgenes = []

for cutoff, topx in zip(cutoffs1, top_x):
    most_prom = above_cutoff(most_common_nodes_in_tis, cutoff)
    most_prom_in_overlap = get_genes_in_overlap(most_prom, interactome_vs_tis_overlap)
    top_x_pathway_hubgenes = get_top_x_genes_in_overlap(pathway_hubgenes, interactome_vs_tis_overlap, topx)
    tis_prom_genes.append(most_prom)
    tis_most_prom_genes_in_int_overlap.append(most_prom_in_overlap)
    tis_top_pathway_hubgenes.append(top_x_pathway_hubgenes)

for cutoff, prom_genes, prom_genes_in_overlap in zip(cutoffs1, tis_prom_genes, tis_most_prom_genes_in_int_overlap):
    print(f"With a cutoff of >= {cutoff}/{len(tis_network_dict)} networks, there are "
          f"{len(prom_genes)} promiscuous genes, "
          f"\n{len(prom_genes_in_overlap)} of which are genes in the overall interactome "
          f"({len(prom_genes_in_overlap) / len(prom_genes):.2%})\n")

fig = plt.figure(figsize=(10, 8))
fig.suptitle('Most connected interactome proteins vs Most common proteins across tissues', weight="bold",
             size="x-large")
i = 1
for pathway_hub, prom_genes_in_overlap, cutoff, topx in zip(tis_top_pathway_hubgenes,
                                                            tis_most_prom_genes_in_int_overlap, cutoffs1, top_x):
    fig.add_subplot(410 + i)
    venn_d(pathway_hub, prom_genes_in_overlap,
           f"Top {topx} most connected interactome proteins",
           f"Most common proteins across tissues w/ cutoff = {cutoff}")
    i += 1
plt.savefig(os.path.join(figures_dir, "most_conn_interactome_vs_most_common_disease.png"))

cutoffs2 = [15, 18, 20, 21]
bottom_x = [100, 200, 300, 400]
tis_least_prom_genes = []
tis_least_prom_genes_in_int_overlap = []
tis_bottom_pathway_hubgenes = []

for cutoff, bottomx in zip(cutoffs2, bottom_x):
    least_prom = below_cutoff(least_common_nodes_in_tis, cutoff)
    least_prom_in_overlap = get_genes_in_overlap(least_prom, interactome_vs_tis_overlap)
    bottom_x_pathway_hubgenes = get_top_x_genes_in_overlap(pathway_hubgenes[::-1], interactome_vs_tis_overlap, bottomx)
    tis_least_prom_genes.append(least_prom)
    tis_least_prom_genes_in_int_overlap.append(least_prom_in_overlap)
    tis_bottom_pathway_hubgenes.append(bottom_x_pathway_hubgenes)

for cutoff, least_prom_genes, least_prom_genes_in_overlap in zip(cutoffs2, tis_least_prom_genes,
                                                                 tis_least_prom_genes_in_int_overlap):
    print(f"With a cutoff of <= {cutoff}/{len(tis_network_dict)} networks, there are "
          f"{len(least_prom_genes)} least promiscuous genes, "
          f"\n{len(least_prom_genes_in_overlap)} of which are genes in the overall interactome "
          f"({len(least_prom_genes_in_overlap) / len(least_prom_genes):.2%})\n")

fig = plt.figure(figsize=(10, 8))
fig.suptitle('Least connected interactome proteins vs Least common proteins across tissues', weight="bold",
             size="x-large")
i = 1
for pathway_hub, l_prom_genes_in_overlap, cutoff, bottomx in zip(tis_bottom_pathway_hubgenes,
                                                                 tis_least_prom_genes_in_int_overlap, cutoffs2,
                                                                 bottom_x):
    fig.add_subplot(410 + i)
    venn_d(pathway_hub, l_prom_genes_in_overlap,
           f"Top {bottomx} least connected interactome proteins",
           f"Least common proteins across tissues w/ cutoff = {cutoff}")
    i += 1
plt.savefig(os.path.join(figures_dir, "least_conn_interactome_vs_least_common_disease.png"))

# Cell types
cutoffs3 = [30, 29, 28, 27]
top_x = [150, 550, 1100, 1800]
ct_prom_genes = []
ct_most_prom_genes_in_int_overlap = []
ct_top_pathway_hubgenes = []

for cutoff, topx in zip(cutoffs3, top_x):
    most_prom = above_cutoff(most_common_nodes_in_ct, cutoff)
    most_prom_in_overlap = get_genes_in_overlap(most_prom, interactome_vs_ct_overlap)
    top_x_pathway_hubgenes = get_top_x_genes_in_overlap(pathway_hubgenes, interactome_vs_ct_overlap, topx)
    ct_prom_genes.append(most_prom)
    ct_most_prom_genes_in_int_overlap.append(most_prom_in_overlap)
    ct_top_pathway_hubgenes.append(top_x_pathway_hubgenes)

for cutoff, prom_genes, prom_genes_in_overlap in zip(cutoffs3, ct_prom_genes, ct_most_prom_genes_in_int_overlap):
    print(f"With a cutoff of >= {cutoff}/{len(ct_network_dict)} networks, there are "
          f"{len(prom_genes)} promiscuous genes, "
          f"\n{len(prom_genes_in_overlap)} of which are genes in the overall interactome "
          f"({len(prom_genes_in_overlap) / len(prom_genes):.2%})\n")

fig = plt.figure(figsize=(10, 8))
fig.suptitle('Most connected interactome proteins vs Most common proteins across cell types', weight="bold",
             size="x-large")
i = 1
for pathway_hub, prom_genes_in_overlap, cutoff, topx in zip(ct_top_pathway_hubgenes, ct_most_prom_genes_in_int_overlap,
                                                            cutoffs3, top_x):
    fig.add_subplot(410 + i)
    venn_d(pathway_hub, prom_genes_in_overlap,
           f"Top {topx} most connected interactome proteins",
           f"Most common proteins across cell types w/ cutoff = {cutoff}")
    i += 1
plt.savefig(os.path.join(figures_dir, "most_conn_interactome_vs_most_common_disease.png"))

cutoffs4 = [9, 10, 12, 13]
bottom_x = [100, 150, 300, 400]
ct_least_prom_genes = []
ct_least_prom_genes_in_int_overlap = []
ct_bottom_pathway_hubgenes = []

for cutoff, bottomx in zip(cutoffs4, bottom_x):
    least_prom = below_cutoff(least_common_nodes_in_ct, cutoff)
    least_prom_in_overlap = get_genes_in_overlap(least_prom, interactome_vs_ct_overlap)
    bottom_x_pathway_hubgenes = get_top_x_genes_in_overlap(pathway_hubgenes[::-1], interactome_vs_ct_overlap, bottomx)
    ct_least_prom_genes.append(least_prom)
    ct_least_prom_genes_in_int_overlap.append(least_prom_in_overlap)
    ct_bottom_pathway_hubgenes.append(bottom_x_pathway_hubgenes)

for cutoff, least_prom_genes, least_prom_genes_in_overlap in zip(cutoffs4, ct_least_prom_genes,
                                                                 ct_least_prom_genes_in_int_overlap):
    print(f"With a cutoff of <= {cutoff}/{len(ct_network_dict)} networks, there are "
          f"{len(least_prom_genes)} least promiscuous genes, "
          f"\n{len(least_prom_genes_in_overlap)} of which are genes in the overall interactome "
          f"({len(least_prom_genes_in_overlap) / len(least_prom_genes):.2%})\n")

fig = plt.figure(figsize=(10, 8))
fig.suptitle('Least connected interactome proteins vs Least common proteins across cell types', weight="bold",
             size="x-large")
i = 1
for pathway_hub, l_prom_genes_in_overlap, cutoff, bottomx in zip(ct_bottom_pathway_hubgenes,
                                                                 ct_least_prom_genes_in_int_overlap, cutoffs4,
                                                                 bottom_x):
    fig.add_subplot(410 + i)
    venn_d(pathway_hub, l_prom_genes_in_overlap,
           f"Top {bottomx} least connected interactome proteins",
           f"Least common proteins across cell types w/ cutoff = {cutoff}")
    i += 1
plt.savefig(os.path.join(figures_dir, "least_conn_interactome_vs_least_common_disease.png"))

# Cell lines
cutoffs5 = [22, 21, 20, 19]
top_x = [100, 500, 1200, 2200]
cl_prom_genes = []
cl_most_prom_genes_in_int_overlap = []
cl_top_pathway_hubgenes = []

for cutoff, topx in zip(cutoffs5, top_x):
    most_prom = above_cutoff(most_common_nodes_in_cl, cutoff)
    most_prom_in_overlap = get_genes_in_overlap(most_prom, interactome_vs_cl_overlap)
    top_x_pathway_hubgenes = get_top_x_genes_in_overlap(pathway_hubgenes, interactome_vs_cl_overlap, topx)
    cl_prom_genes.append(most_prom)
    cl_most_prom_genes_in_int_overlap.append(most_prom_in_overlap)
    cl_top_pathway_hubgenes.append(top_x_pathway_hubgenes)

for cutoff, prom_genes, prom_genes_in_overlap in zip(cutoffs5, cl_prom_genes, cl_most_prom_genes_in_int_overlap):
    print(f"With a cutoff of >= {cutoff}/{len(cl_network_dict)} networks, there are "
          f"{len(prom_genes)} promiscuous genes, "
          f"\n{len(prom_genes_in_overlap)} of which are genes in the overall interactome "
          f"({len(prom_genes_in_overlap) / len(prom_genes):.2%})\n")

fig = plt.figure(figsize=(10, 8))
fig.suptitle('Most connected interactome proteins vs Most common proteins across cell lines', weight="bold",
             size="x-large")
i = 1
for pathway_hub, prom_genes_in_overlap, cutoff, topx in zip(cl_top_pathway_hubgenes, cl_most_prom_genes_in_int_overlap,
                                                            cutoffs5, top_x):
    fig.add_subplot(410 + i)
    venn_d(pathway_hub, prom_genes_in_overlap,
           f"Top {topx} most connected interactome proteins",
           f"Most common proteins across cell lines w/ cutoff = {cutoff}")
    i += 1
plt.savefig(os.path.join(figures_dir, "most_conn_interactome_vs_most_common_disease.png"))

cutoffs6 = [9, 11, 12, 13]
bottom_x = [100, 200, 300, 400]
cl_least_prom_genes = []
cl_least_prom_genes_in_int_overlap = []
cl_bottom_pathway_hubgenes = []

for cutoff, bottomx in zip(cutoffs6, bottom_x):
    least_prom = below_cutoff(least_common_nodes_in_cl, cutoff)
    least_prom_in_overlap = get_genes_in_overlap(least_prom, interactome_vs_cl_overlap)
    bottom_x_pathway_hubgenes = get_top_x_genes_in_overlap(pathway_hubgenes[::-1], interactome_vs_cl_overlap, bottomx)
    cl_least_prom_genes.append(least_prom)
    cl_least_prom_genes_in_int_overlap.append(least_prom_in_overlap)
    cl_bottom_pathway_hubgenes.append(bottom_x_pathway_hubgenes)

for cutoff, least_prom_genes, least_prom_genes_in_overlap in zip(cutoffs6, cl_least_prom_genes,
                                                                 cl_least_prom_genes_in_int_overlap):
    print(f"With a cutoff of <= {cutoff}/{len(cl_network_dict)} networks, there are "
          f"{len(least_prom_genes)} least promiscuous genes, "
          f"\n{len(least_prom_genes_in_overlap)} of which are genes in the overall interactome "
          f"({len(least_prom_genes_in_overlap) / len(least_prom_genes):.2%})\n")

fig = plt.figure(figsize=(10, 8))
fig.suptitle('Least connected interactome proteins vs Least common proteins across cell types', weight="bold",
             size="x-large")
i = 1
for pathway_hub, l_prom_genes_in_overlap, cutoff, bottomx in zip(cl_bottom_pathway_hubgenes,
                                                                 cl_least_prom_genes_in_int_overlap, cutoffs6,
                                                                 bottom_x):
    fig.add_subplot(410 + i)
    venn_d(pathway_hub, l_prom_genes_in_overlap,
           f"Top {bottomx} least connected interactome proteins",
           f"Least common proteins across cell types w/ cutoff = {cutoff}")
    i += 1
plt.savefig(os.path.join(figures_dir, "least_conn_interactome_vs_least_common_disease.png"))
