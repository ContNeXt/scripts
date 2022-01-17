# -*- coding: utf-8 -*-

"""Pathway-level analysis for ContNeXt"""

import json
import os
from itertools import product

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from clustergrammer2 import Network
from tqdm import tqdm

import ora
from network_utils import create_network_from_edge_file, edge_file_path, pathway_similarity, df_and_normalize_dict

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
                   desc="Creating/loading tissue network objects") if ID != ".DS_Store"
}

# cell type networks
ct_network_dict = {
    ID: create_network_from_edge_file(edge_file_path(os.path.join(data_dir, "coexpr_networks", "cell_type"), ID),
                                      ID)
    for ID in tqdm(os.listdir(os.path.join(data_dir, "coexpr_networks", "cell_type")),
                   desc="Creating/loading cell type network objects") if ID != ".DS_Store"
}

# cell line  networks
cl_network_dict = {
    ID: create_network_from_edge_file(edge_file_path(os.path.join(data_dir, "coexpr_networks", "cell_line"), ID),
                                      ID)
    for ID in tqdm(os.listdir(os.path.join(data_dir, "coexpr_networks", "cell_line")),
                   desc="Creating/loading cell line network objects") if ID != ".DS_Store"
}

# KEGG pathway assignments
gene_pathway_file = os.path.join(data_dir, "pathway", "gene_pathway_assignment.json")
with open(gene_pathway_file, 'r') as f:
    gene_pathways_dict = json.load(f)

pathway_genes = os.path.join(data_dir, "pathway", "kegg_hgnc_ids.gmt")
pathway_genes_dict = ora.gmt_parser(pathway_genes, 3, 10000)

# ontology name mappings
with open(os.path.join(data_dir, "mappings", "uberon_name_mappings.json"), 'r') as f:
    uberon_name_mappings = json.load(f)
with open(os.path.join(data_dir, "mappings", "CL_name_mappings.json"), 'r') as f:
    CL_name_mappings = json.load(f)
with open(os.path.join(data_dir, "mappings", "CLO_name_mappings.json"), 'r') as f:
    CLO_name_mappings = json.load(f)

# pathway names mapping
with open(os.path.join(data_dir, "pathway", 'kegg_mapping.json'), 'r') as f:
    as_str = f.read()
    pathway_name_mapping = eval(as_str)

pathway_name_mapping = {path.split(":")[1]: name.split(" - ")[0] for path, name in pathway_name_mapping.items()}
pathway_id_mapping = {v: k for k, v in pathway_name_mapping.items()}

# ### Create matrices and clustergrammer objects ####

blue_color_palette = sns.color_palette("Blues", as_cmap=True)


def plot_df_with_clustering(df, context_name, figname=None):
    mpl.rcParams['xtick.labelsize'] = 'x-small'
    mpl.rcParams['ytick.labelsize'] = 'small'
    g = sns.clustermap(df, cmap=blue_color_palette, figsize=(40, 10),  # row_cluster=False,
                       cbar_pos=None,
                       yticklabels=df.index,
                       xticklabels=df.columns)
    g.ax_col_dendrogram.set_visible(False)
    g.ax_row_dendrogram.set_visible(False)
    plt.title(f"Pathway fingerprints per {context_name}", fontsize=45, weight="bold")
    plt.xlabel("Pathway name", fontsize=30)
    plt.ylabel(f"{context_name.capitalize()} name", fontsize=30)
    if figname:
        plt.savefig(figname, bbox_inches='tight', dpi=320)
    return g


# Tissues
pathway_similarity_per_tis = {ID: {pathway: None for pathway in pathway_genes_dict.keys()}
                              for ID in tis_network_dict.keys()}

for ID, (pathway, gene_list) in tqdm(product(tis_network_dict.keys(), pathway_genes_dict.items()),
                                     desc="Testing each pathway on each tissue network",
                                     total=len(pathway_genes_dict) * len(tis_network_dict)):
    is_enriched = pathway_similarity(gene_list, tis_network_dict[ID])
    pathway_similarity_per_tis[ID][pathway] = is_enriched

tissue_df = df_and_normalize_dict(pathway_similarity_per_tis)
tissue_df.rename(index=uberon_name_mappings, columns=pathway_name_mapping, inplace=True)

tissue_heatmap = plot_df_with_clustering(tissue_df, "tissue", os.path.join(figures_dir, "tissue_vs_pathways_similarity.pdf"))

tis_net = Network()
tis_net.load_df(tissue_df)
tis_net.cluster()
# save visualization JSON to file for use by front end
tis_net.write_json_to_file('viz', os.path.join(data_dir, "misc_data", 'tissue-pathway_clustergrammer.json'))

# Cell types
pathway_similarity_per_ct = {ID: {pathway: None for pathway in pathway_genes_dict.keys()}
                             for ID in ct_network_dict.keys()}

for ID, (pathway, gene_list) in tqdm(product(ct_network_dict.keys(), pathway_genes_dict.items()),
                                     desc="Testing each pathway on each cell type network",
                                     total=len(pathway_genes_dict) * len(ct_network_dict)):
    is_enriched = pathway_similarity(gene_list, ct_network_dict[ID])
    pathway_similarity_per_ct[ID][pathway] = is_enriched

celltype_df = df_and_normalize_dict(pathway_similarity_per_ct)
celltype_df.rename(index=CL_name_mappings, columns=pathway_name_mapping, inplace=True)

celltype_heatmap = plot_df_with_clustering(celltype_df, "cell type", os.path.join(figures_dir, "celltype_vs_pathways_similarity.pdf"))

ct_net = Network()
ct_net.load_df(celltype_df)
ct_net.cluster(dist_type="euclidean")
# save visualization JSON to file for use by front end
ct_net.write_json_to_file('viz', os.path.join(data_dir, "misc_data", 'celltype-pathway_clustergrammer.json'))

# cell lines
pathway_similarity_per_cl = {ID: {pathway: None for pathway in pathway_genes_dict.keys()}
                             for ID in cl_network_dict.keys()}

for ID, (pathway, gene_list) in tqdm(product(cl_network_dict.keys(), pathway_genes_dict.items()),
                                     desc="Testing each pathway on each cell line network",
                                     total=len(pathway_genes_dict) * len(cl_network_dict)):
    is_enriched = pathway_similarity(gene_list, cl_network_dict[ID])
    pathway_similarity_per_cl[ID][pathway] = is_enriched

cellline_df = df_and_normalize_dict(pathway_similarity_per_cl)
cellline_df.rename(index=CLO_name_mappings, columns=pathway_name_mapping, inplace=True)

cellline_heatmap = plot_df_with_clustering(cellline_df, "cell line", os.path.join(figures_dir, "cellline_vs_pathways_similarity.pdf"))

cl_net = Network()
cl_net.load_df(cellline_df)
cl_net.cluster(dist_type="euclidean")
# save visualization JSON to file for use by front end
cl_net.write_json_to_file('viz', os.path.join(data_dir, "misc_data", 'cellline-pathway_clustergrammer.json'))
