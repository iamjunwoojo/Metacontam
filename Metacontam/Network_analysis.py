import subprocess
import os
import pandas as pd
import numpy as np
import networkx as nx
import pickle
import matplotlib.pyplot as plt
from networkx.algorithms import community
from Metacontam.louvain.louvain import detect_communities
from Metacontam.louvain.modularity import modularity
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to save figures to file



def Make_adjacent_matrix(Rscript_path , input_file ,output, coefficient_threshold):
    input_file = os.path.join(output,"kraken_filtered_matrix.txt")
    output_file = os.path.join(output,"Network_Output","edge.tsv")
    threshold = str(coefficient_threshold)
    subprocess.run(
        ["Rscript", Rscript_path, input_file, output_file,threshold],
        check=True
    )


def community_detection(output):
    edge_df=pd.read_csv(os.path.join(output,"Network_Output","edge.tsv"),sep="\t")
    node_input=np.unique(np.array(edge_df['v1'].to_list()+edge_df['v2'].to_list(),dtype=str))
    edge_list=[]
    # for line in edge_df.iterrows():
    #     v1=str(int(line[1]['v1']))
    #     v2=str(int(line[1]['v2']))
    #     asso=float(line[1]['asso'])
    #     edge_list.append((v1,v2,asso))
    for line in edge_df.iterrows():
        v1=str(int(line[1]['v1']))
        v2=str(int(line[1]['v2']))
        asso=float(line[1]['asso'])
        if asso >0:
            edge_list.append((v1,v2,asso))


    edge_list=[(int(i),int(l),z) for i,l,z in edge_list]
    node_input=[int(i) for i in node_input]


    G = nx.Graph()
    G.add_nodes_from(node_input)
    G.add_weighted_edges_from(edge_list)
    return G




def draw_small_communities(G, high_preval, output, node_size=20, alpha=1, k=None, randomized=False, seed=5):
    # Detect communities in the graph
    partition = detect_communities(G, high_preval, randomized=randomized, verbose=False)
    print("Modularity for best partition:", modularity(G, partition))

    # Save the partition to a file
    partition_path = os.path.join(output, "Network_Output", "partition.pkl")
    os.makedirs(os.path.dirname(partition_path), exist_ok=True)  # Create directories if they don't exist
    with open(partition_path, "wb") as f:
        pickle.dump(partition, f)
    print(f"Partition saved to {partition_path}")

    # Select only communities with more than 10 nodes
    community_map = {}
    for community, nodes in enumerate(partition):
        if len(nodes) > 10:  # Include communities with more than 10 nodes
            for node in nodes:
                community_map[node] = community

    # Filter nodes to plot only those in the selected communities
    small_community_nodes = set(community_map.keys())
    small_community_subgraph = G.subgraph(small_community_nodes)

    # Generate the graph layout and define color mapping
    cmap = plt.get_cmap("Dark2")
    pos = nx.spring_layout(small_community_subgraph, k=k, seed=seed)
    indexed = [community_map.get(node) for node in small_community_subgraph]

    # Plot the graph
    plt.figure(figsize=(10, 8))
    plt.axis("off")
    nx.draw_networkx_nodes(
        small_community_subgraph, pos=pos, cmap=cmap, node_color=indexed,
        node_size=node_size, alpha=alpha
    )
    nx.draw_networkx_edges(small_community_subgraph, pos=pos, alpha=0.2)

    # Add labels to the nodes
    labels = {node: str(node) for node in small_community_subgraph.nodes()}
    nx.draw_networkx_labels(small_community_subgraph, pos=pos, labels=labels, font_size=8)

    # Save the graph to the specified path
    output_path = os.path.join(output, "Network_Output", "Figure")
    os.makedirs(output_path, exist_ok=True)  # Create directories if they don't exist
    plt.savefig(os.path.join(output_path, "Community.png"), format="png", dpi=300)
    plt.close()

def contam_community(partition, high_preval):
    for community in partition:
        if set(community).intersection(high_preval)==set(high_preval):
            return community
