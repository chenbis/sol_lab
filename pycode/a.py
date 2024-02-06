import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# Read the distance matrix from CSV file
distance_matrix_df = pd.read_csv("distance_matrix.csv", index_col=0)

# Example amino acid sequence (replace this with your actual sequence)
sequences = pd.read_csv('files/forchen_F_26L.csv')["cdr3_amino_acid"].tolist()

neighbors = {}


def calculate_distance(seq1, seq2):
    distance = 0
    for aa1, aa2 in zip(seq1, seq2):
        # Assuming aa1 and aa2 are indices for the amino acids
        distance += distance_matrix_df[aa1][aa2]
    return distance

# Step 3: Find Closest Sequence
for i, seq1 in enumerate(sequences):
    min_distance = float('inf')
    closest_sequence_index = -1
    
    for j, seq2 in enumerate(sequences):
        if i != j:  # Avoid comparing a sequence with itself
            current_distance = calculate_distance(seq1, seq2)
            
            if current_distance < min_distance:
                min_distance = current_distance
                closest_sequence_index = j
    
    closest_sequence = sequences[closest_sequence_index]
    
    # print(f"Sequence {sequences[i]}: Closest Sequence is {sequences[closest_sequence_index]} with distance {min_distance}")
    neighbors[sequences[i]] = sequences[closest_sequence_index]


# Create a graph
G = nx.Graph()

# Add nodes and edges to the graph
for node, neighbors in neighbors.items():
    G.add_node(node)
    G.add_edges_from((node, neighbor) for neighbor in neighbors)

# Draw the graph
pos = nx.spring_layout(G)  # You can use other layout algorithms as well
nx.draw(G, pos, with_labels=True, font_weight='bold', node_size=700, node_color='skyblue', font_size=10, edge_color='gray', width=1.5)

# Show the graph
plt.show()