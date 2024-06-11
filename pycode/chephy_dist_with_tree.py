import pandas as pd
from collections import defaultdict
import csv
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

def truncate_sequences(csv_file, column_name, right, left):

    sequences_set = set()
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sequence = row[column_name].strip()
            if sequence == "":
                continue
            start_index = max(0, len(sequence)//2 - left)
            end_index = min(len(sequence), len(sequence)//2 + right)
            sequences_set.add(sequence[start_index:end_index])
            
    return sequences_set

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

def find_che_phy_dist(sequences_set, max_mutations):
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

    for seq in neighbors:
        for var in neighbors[seq]:
            couples[seq].append([var, sum(distances_dict[(aa1, aa2)] for aa1, aa2 in zip(seq, var))])

    couples = {key: sorted(value, key=lambda x: x[1]) for key, value in couples.items()}
    return couples




def write_couples_file(couples, directory, filename):
    path = Path(directory)

    path.mkdir(parents=True, exist_ok=True)
    output_file = "{}/{}.json".format(directory, filename)
    with open(output_file, "w") as f:
        ujson.dump(couples, f)

def main():

    max_mutations = 3
    out_folder = "output_files/tests"
    out_name = "test"

    # 1000 sequences
    sequences_file = 'files/forchen_F_26L.csv'
    sequences_headers = ["cdr3_amino_acid"]

    # # 40K sequences
    # sequences_file = 'files/for_chen_B2.csv'
    # sequences_headers = ["CDR3.aa"]

    # # 1M sequences
    # sequences_file = "random_sequences.csv"
    # sequences_headers = ["sequences"]

    # sequences_file = 'files/articles/neo.csv'
    # sequences_headers = ["Melan.A", "MAGE.A3", "GP100","MAGE.A10", "PHLPP2", "HHAT", "ZCCHC11",\
    #                      "SCL25A48", "MUC1", "MMP9", "UTP20", "Influenza.A.PB1",\
    #                         "EBV.BMLF1","hCMV.pp65","hCMV.pp65.2"]


    for header in sequences_headers:
        sequences_set = truncate_sequences(sequences_file, header, 4, 4)
        couples = find_che_phy_dist(sequences_set, max_mutations)
        write_couples_file(couples, "{}/{}".format(out_folder, header), out_name)

if __name__ == "__main__":
    main()