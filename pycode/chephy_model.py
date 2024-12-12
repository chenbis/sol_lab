import pandas as pd
from collections import defaultdict
from tqdm import tqdm
from pathlib import Path
import ujson

def create_tree(sequences):
    tree = {}
    for sequence in sequences:
        current_dict = tree
        for letter in sequence:
            current_dict = current_dict.setdefault(letter, {})
    return tree

def traverse_tree(tree, sequence, max_diff, path="", current_diff=0, depth=0):
    if current_diff > max_diff:
        return []
    
    if not tree:
        return [path] if current_diff <= max_diff else []
    
    matching_sequences = []
    for char, subtree in tree.items():
        new_diff = current_diff + (1 if (depth < len(sequence) and char != sequence[depth]) else 0)
        matching_sequences.extend(traverse_tree(subtree, sequence, max_diff, path + char, new_diff, depth + 1))

    return matching_sequences

def find_sequences_within_distance(tree, sequence, max_diff):
    return traverse_tree(tree, sequence, max_diff)

def truncate_sequences(sequences, right=4, left=4):

    full_to_trunc_map = defaultdict(set)

    sequences_set = set()
    for sequence in sequences:
        if sequence == "":
            continue
        start_index = max(0, len(sequence)//2 - left)
        end_index = min(len(sequence), len(sequence)//2 + right)
        trunc_seq = sequence[start_index:end_index]
        sequences_set.add(trunc_seq)
        full_to_trunc_map[trunc_seq].add(sequence)

    return sequences_set, full_to_trunc_map

def invert_dict(d): 
    inverse = dict() 
    for key in d: 
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse: 
                # If not create a new list
                inverse[item] = {key} 
            else: 
                inverse[item].add(key) 
    return inverse

def find_che_phy_dist(sequences_set, max_mutations, max_dist):
    """
    finds chemical physical distance between all sequences in the set
    sequences set - set of sequences to compare
    max_mutations - maximum number of different amino acids between two close sequences
    """

    neighbors = {}
    tree = create_tree(sequences_set)
    for seq in tqdm(sequences_set):
        results = find_sequences_within_distance(tree, seq, max_mutations)
        results.remove(seq)
        if results:
            neighbors[seq] = set(results)
    
    
    inverted = invert_dict(neighbors)
    neighbors = defaultdict(set, {k: neighbors.get(k, set()) | inverted.get(k, set())\
                                   for k in set(neighbors) | set(inverted)})
    

    distances_csv = "distance_matrix.csv"
    distances_df = pd.read_csv(distances_csv, index_col=0)
    distances_dict = {(aa1, aa2): distances_df.loc[aa1, aa2] for aa1 in distances_df.index for aa2 in distances_df.columns}
    
    couples = defaultdict(list)
    worst_case_distances = {seq: get_worst_case_distance(seq, distances_df) for seq in neighbors}
    # Optimize the loop
    for seq in tqdm(neighbors):
        seq_neighbors = neighbors[seq]  # Avoid repeated lookup
        worst_case_distance = worst_case_distances[seq]  # Retrieve precomputed value

        for var in seq_neighbors:
            # Sum distances using zip
            actual_distance = sum(distances_dict[(aa1, aa2)] for aa1, aa2 in zip(seq, var))
            
            # Normalize the distance
            normalized_distance = actual_distance / worst_case_distance if worst_case_distance > 0 else 0
            normalized_distance = float('%.3f'%(normalized_distance))
            if not normalized_distance > max_dist:
                couples[seq].append([var, normalized_distance])

    couples = {key: sorted(value, key=lambda x: x[1]) for key, value in couples.items()}
    return couples

def write_couples_file(couples, directory, filename):
    path = Path(directory)

    path.mkdir(parents=True, exist_ok=True)
    output_file = "{}/{}.json".format(directory, filename)
    with open(output_file, "w") as f:
        ujson.dump(couples, f)



def get_worst_case_distance(seq, substitution_matrix):
    max_distance = 0
    for aa in seq:
        # For each amino acid, find the maximum possible substitution cost
        max_substitution_cost = substitution_matrix[aa].max()
        max_distance += max_substitution_cost
    return max_distance