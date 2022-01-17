# -*- coding: utf-8 -*-

"""Network overview analysis for ContNeXt"""

import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

import ora
from network_utils import create_network_from_edge_file, edge_file_path, load_interactome, most_common

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
    ID: create_network_from_edge_file(edge_file_path(os.path.join(data_dir, "coexpr_networks", "cell.line"), ID),
                                      ID)
    for ID in tqdm(os.listdir(os.path.join(data_dir, "coexpr_networks", "cell.line")),
                   desc="Creating/loading network objects") if ID != ".DS_Store"
}

# interactome PPI network
interactome = load_interactome(os.path.join(data_dir, "interactome", "interactome_18_01_2021.tsv"))

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
with open(os.path.join(data_dir, "mappings", "hgnc_name_mappings.json"), 'r') as f:
    hgnc_mappings = json.load(f)

# group context data into one variable¶
contexts = {
    "tissues": {
        "nets": tis_network_dict,
        "name_mapping": uberon_name_mappings
    },
    "cell types": {
        "nets": ct_network_dict,
        "name_mapping": CL_name_mappings
    },
    "cell lines": {
        "nets": cl_network_dict,
        "name_mapping": CLO_name_mappings
    }
}
IDs_per_context = {
    group: sorted(group_data["nets"].keys(), key=lambda ID: len(group_data["nets"][ID].nodes()), reverse=True)
    for group, group_data in contexts.items()}

# ### LOAD/PROCESS CONTROLLABILITY ANALYSIS RESULTS ####
summary_keys = [
    "Number of nodes(N)",
    "Number of edges(E)",
    "Average degree(Kmean)",
    "Number of driver nodes(Nd)",
    "Fraction of driver nodes(Nd/N)",
    "Fraction of type-I critical nodes",
    "Fraction of type-I redundant nodes",
    "Fraction of type-I ordinary nodes",
    "Fraction of type-II critical nodes",
    "Fraction of type-II redundant nodes",
    "Fraction of type-II ordinary nodes",
    "Fraction of critical links",
    "Fraction of redundant links",
    "Fraction of ordinary links",
    "Average local clustering coefficient",
    "Global clustering coefficient"
    "Average directed local clustering coefficient",
    "Global directed clustering coefficient",
    "Average betweenness centralities",
    "Average closeness centralities",
    "Average authority centralities",
    "Average hub centralities",
    "Fraction of source nodes",
    "Fraction of external dilations",
    "Fraction of internal dilations",
]

with open(os.path.join(data_dir, "controllability_analysis", "interactome.output"), "r") as f:
    summary_raw = f.read().strip().split(",")

summary = dict(zip(summary_keys, summary_raw))
with open(os.path.join(data_dir, "controllability_analysis", "interactome_summary.tsv"), 'w') as f:
    for key in summary.keys():
        f.write(f"{key}\t{summary[key]}\n")

node_type = pd.read_table(os.path.join(data_dir, "controllability_analysis", "interactome.nodetype"), sep=" ", usecols=["#Name", "TypeI"])
node_type['#Name'] = node_type['#Name'].astype(str)
importance = {0: "critical", 1: "redundant", 2: "ordinary"}
node_type = node_type.replace({"#Name": hgnc_mappings, "TypeI": importance})
node_type.columns = ["protein", "classification"]
node_type.to_csv(os.path.join(data_dir, "controllability_analysis", "interactome_node_classifications.tsv"), sep='\t', index=False)

indispensable_nodes = list(node_type.loc[node_type.classification == "critical"][
                               "protein"])  # increases driver nodes if removed ex. in a directed path
dispensable_nodes = list(node_type.loc[node_type.classification == "redundant"][
                             "protein"])  # decreases driver nodes if removed ex. leaf in a star
neutral_nodes = list(node_type.loc[node_type.classification == "ordinary"][
                         "protein"])  # no change in driver nodes if removed ex. central hub of a star
print(f"there are {len(indispensable_nodes)} indispensable nodes, {len(dispensable_nodes)} dispensable nodes, and {len(neutral_nodes)} neutral nodes.")

with open(os.path.join(data_dir, "controllability_analysis", "interactome_indispensable_nodes.txt"), 'w') as f:
    f.write("\n".join(indispensable_nodes))

edge_type = pd.read_table(os.path.join(data_dir, "controllability_analysis", "interactome.linktype"), sep=" ")
edge_type['#source'] = edge_type['#source'].astype(str)
edge_type['target'] = edge_type['target'].astype(str)
edge_type = edge_type.replace({"#source": hgnc_mappings, "target": hgnc_mappings, "LinkType": importance})
edge_type.columns = ["source", "target", "classification"]
edge_type.to_csv(os.path.join(data_dir, "controllability_analysis", "interactome_edge_classifications.tsv"), sep='\t', index=False)

indispensable_edges = list(zip(edge_type.loc[edge_type.classification == "critical"]["source"],
                               edge_type.loc[edge_type.classification == "critical"][
                                   "target"]))  # increases driver nodes if removed
dispensable_edges = list(zip(edge_type.loc[edge_type.classification == "redundant"]["source"],
                             edge_type.loc[edge_type.classification == "redundant"][
                                 "target"]))  # decreases driver nodes if removed
neutral_edges = list(zip(edge_type.loc[edge_type.classification == "ordinary"]["source"],
                         edge_type.loc[edge_type.classification == "ordinary"][
                             "target"]))  # no change in driver nodes if removed
print(f"there are {len(indispensable_edges)} indispensable edges, {len(dispensable_edges)} dispensable edges, and {len(neutral_edges)} neutral edges.")

pd.DataFrame(indispensable_edges).to_csv(os.path.join(data_dir, "controllability_analysis", "interactome_indispensable_edges.txt"),
                                         sep='\t', index=False, header=False)

most_common_nodes_in_tis = most_common(tis_network_dict.values(), comparison="nodes")
most_common_nodes_in_ct = most_common(ct_network_dict.values(), comparison="nodes")
most_common_nodes_in_cl = most_common(cl_network_dict.values(), comparison="nodes")

# ### FIGURES ####

# Network size overview figure
term_and_node_num_per_context = {}
# use a count to preserve ordering
count = 0
for group, ids in IDs_per_context.items():
    for ID in ids:
        term_and_node_num_per_context[count] = {
            "Context": group,
            "Term": contexts[group]["name_mapping"][ID].capitalize(),
            "Number of nodes": len(contexts[group]["nets"][ID].nodes()),
        }
        count += 1
node_count_df = pd.DataFrame(term_and_node_num_per_context.values(), index=term_and_node_num_per_context.keys())

sns.reset_defaults()
f, ax = plt.subplots(figsize=(20, 15))
# the node_count_df has the preserved intended order
legend_order = node_count_df["Context"].unique()
sns.barplot(x="Number of nodes", y="Term", hue="Context", hue_order=legend_order, data=node_count_df, dodge=False)
# extend chart above and below
ax.set_ylim([len(node_count_df), -1.5])
plt.title("Network size distribution", fontsize=45, x=0.25)
ax.set_xlabel(ax.get_xlabel(), fontsize=25)
ax.set_ylabel("Individual class and corresponding contexts", fontsize=25)
for p in ax.patches:
    # label bars
    ax.annotate(int(np.nan_to_num(p.get_width(), 0)), xy=(p.get_width(), p.get_y() + p.get_height() / 2),
                xytext=(5, 0), textcoords='offset points', ha="left", va="center")
# move legend to outside
plt.legend(bbox_to_anchor=(1.01, 0.9), loc=2, borderaxespad=0., fontsize='large')
# save figure
plt.savefig(os.path.join(figures_dir, "Network_Size_Distribution.png"), bbox_inches='tight', dpi=320)

# Term overview figure
term_overview = {
    group: {
        "ids": list(group_data["nets"].keys()),
        "terms": [group_data["name_mapping"][ID] for ID in group_data["nets"].keys()]
    }
    for group, group_data in contexts.items()
}
naming_for_files = {"tissues": "tissue", "cell types": "celltype", "cell lines": "cellline"}
context_info = {}
for context in term_overview:
    with open(os.path.join(data_dir, "misc_data", f"FULL_{naming_for_files[context]}_overview_after_download.tsv"),
              "r") as f:
        lines = [f.readline().strip()]
        lines += [line.strip() for line in f.readlines() if
                  line.strip().split("\t")[1] in term_overview[context]["terms"]]
    context_info[context] = {
        line.split("\t")[0].split(":")[1]: {
            "Name": line.split("\t")[1],
            "Number of datasets": int(line.split("\t")[2]),
            "Number of samples": int(line.split("\t")[3]),
            "Number of nodes": len(contexts[context]["nets"][line.split("\t")[0].split(":")[1]].nodes()),
            "Number of edges": len(contexts[context]["nets"][line.split("\t")[0].split(":")[1]].edges()),
        }
        for line in lines[1:]
    }
    with open(os.path.join(data_dir, "misc_data", f"{naming_for_files[context]}_overview.tsv"), "w") as f:
        f.write("\n".join(lines))

# Correlation between number of nodes to number of datasets and number of samples for all terms figure
tissue_df = pd.DataFrame.from_dict(context_info["tissues"], orient="index")
tissue_df["Context"] = ["tissue"] * len(tissue_df)
celltype_df = pd.DataFrame.from_dict(context_info["cell types"], orient="index")
celltype_df["Context"] = ["cell type"] * len(celltype_df)
cellline_df = pd.DataFrame.from_dict(context_info["cell lines"], orient="index")
cellline_df["Context"] = ["cell line"] * len(cellline_df)
term_info_df = pd.concat([tissue_df, celltype_df, cellline_df], axis=0)
mpl.rcParams['xtick.labelsize'] = 'xx-small'
mpl.rcParams['ytick.labelsize'] = 'xx-small'
ax = sns.pairplot(term_info_df,
                  x_vars=["Number of datasets", "Number of samples"], y_vars=["Number of nodes"],
                  hue="Context", plot_kws={"s": 15, 'alpha': 0.8})
plt.title("Correlation between number of datasets/samples and network size", x=0)
ax._legend.remove()
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., prop={'size': 8})
plt.savefig(os.path.join(figures_dir, "datasets-samples_vs_nodes.png"), bbox_inches='tight', dpi=320)

# Distribution of the frequency of all proteins across the co-expression networks
fig, axes = plt.subplots(1, 3, sharey=True, figsize=(20, 5))
gene_counts_tis = [int(count[0].split('/')[0]) for gene, count in most_common_nodes_in_tis]
gene_counts_ct = [int(count[0].split('/')[0]) for gene, count in most_common_nodes_in_ct]
gene_counts_cl = [int(count[0].split('/')[0]) for gene, count in most_common_nodes_in_cl]

sns.histplot(gene_counts_tis, kde=True, bins=26, ax=axes[0])
sns.histplot(gene_counts_ct, kde=True, bins=26, ax=axes[1])
sns.histplot(gene_counts_cl, kde=True, bins=26, ax=axes[2])

axes[0].set(xlabel="tissues", ylabel="")
axes[1].set(xlabel="cell types", ylabel="")
axes[2].set(xlabel="cell lines", ylabel="")

fig.text(.12, .95, "Protein frequency distribution in co-expression networks", fontsize=35, weight="bold")
fig.text(0.5, -0.03, 'Number of networks', ha='center', fontsize=20)
fig.text(0.08, 0.5, 'Protein count', va='center', rotation='vertical', fontsize=20)

plt.savefig(os.path.join(figures_dir, "prot-freq-dist.png"), bbox_inches='tight', dpi=320)
