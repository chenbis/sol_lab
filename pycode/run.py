import chephy_model as cpm
from collections import defaultdict
import time
import networkx as nx
import pandas as pd
from tqdm import tqdm
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score

def map_trunc_to_full(couples, full_to_trunc_map):
    # Initialize an empty dictionary for the full couples
    couples_full = defaultdict(list)

    # Iterate over each key, value pair in the couples dictionary
    for short_seq, neighbors in couples.items():
        # Get the original sequences for the current short sequence
        original_seqs = set(full_to_trunc_map.get(short_seq, []))

        # Iterate over each original sequence
        for original_seq in original_seqs:
            
            # Iterate over each neighbor in the neighbors list
            for neighbor, distance in neighbors:
                # Get the original sequences for the neighbor
                neighbor_original_seqs = full_to_trunc_map.get(neighbor, [])
                
                # Add each original neighbor sequence with the same distance
                for neighbor_original_seq in neighbor_original_seqs:
                    couples_full[original_seq].append([neighbor_original_seq, distance])

    return couples_full


def find_close_sequences(cdr3, max_mutations=3, right=4, left=4):
    sequences_set, full_to_trunc_map = cpm.truncate_sequences(cdr3, 4, 4)
    couples_trunc = cpm.find_che_phy_dist(sequences_set, max_mutations)
    couples_full = map_trunc_to_full(couples_trunc, full_to_trunc_map)
    return couples_full


def map_clusters(couples, max_neig):
    # Create a graph
    G = nx.Graph()

    # Add edges to the graph based on the dictionary
    for node, neighbors in couples.items():
        count = 0
        for neighbor, weight in neighbors:
            if count < max_neig:
                G.add_edge(node, neighbor, weight=weight)
                count += 1


    # Find connected components
    clusters = list(nx.connected_components(G))

    # Create a dictionary to map sequences to cluster numbers
    cluster_dict = {}
    # for cluster_number, cluster in enumerate(clusters):
    #     for sequence in cluster:
    #         cluster_dict[sequence] = cluster_number    

    for cluster_number, cluster in enumerate(clusters):
        cluster_dict[cluster_number] = cluster
   
    return cluster_dict



def get_true_clusters(data, cdr3_header, eptiope_header):
    clusters_dict = {}

    # Assign clusters based on common antigen
    for cluster_id, (antigen, group) in enumerate(data.groupby(eptiope_header)):
        for cdr3 in group[cdr3_header]:
            clusters_dict[cdr3] = cluster_id

    return clusters_dict

def classify(clusters, data):
    cluster_classification = {}
    for cluster, tcrs in clusters.items():
        epitope_strengths = defaultdict(int)
        for tcr in list(tcrs):
            if tcr == "CASSEGWHSYEQYF":
                print()
            row = data[data['cdr3'] == tcr]
            epitope = row['antigen.epitope'].values[0]
            strength = row['vdjdb.score'].values[0]
            epitope_strengths[epitope] += strength

        best_epitope = max(epitope_strengths, key=epitope_strengths.get)
        cluster_classification[cluster] = best_epitope

    return cluster_classification




def hamming_distance(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def find_sequences_within_distance(cdr3_list, max_dist):
    """Find all sequences that are within a Hamming distance of 1 for each sequence."""
    result = defaultdict(list)
    
    for i, seq1 in enumerate(cdr3_list):
        # Initialize an empty list for each sequence
        
        for seq2 in cdr3_list:
            # Skip comparing the sequence with itself
            if seq1 == seq2:
                continue
            
            # If the Hamming distance is 1, add to the list
            distance = hamming_distance(seq1, seq2)
            if distance <= max_dist:
                result[seq1].append([seq2, distance/8])
    
    return result


def main():

    max_mutations = 8
    out_folder = "output_files/tests"
    out_name = "test"
    max_neig = 1

    # # 1000 sequences
    # data = pd.read_csv('files/forchen_F_26L.csv')
    # cdr3_header = "cdr3_amino_acid"


    # vdjdb beta chain
    data = pd.read_csv("files/vdjdb_cdr3.csv")
    data = data[data["vdjdb.score"] == 3]
    
    cdr3_header = "cdr3"
    epitope_header = "antigen.epitope"

    data.drop_duplicates(subset=cdr3_header)
    cdr3 = list(data[cdr3_header])
    
    couples_full = find_close_sequences(cdr3, max_mutations)
    
    start_time = time.time()
    cpm.write_couples_file(couples_full, "{}/{}".format(out_folder, cdr3_header), "{}_full".format(out_name))
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print(f"Time taken to write full couples: {elapsed_time:.2f} minutes")


    start_time = time.time()
    cluster_dict_pred = map_clusters(couples_full, max_neig)
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print(f"{elapsed_time:.2f} minutes")
    # cluster_dict_true = get_true_clusters(data, cdr3_header, epitope_header)

    start_time = time.time()
    cluster_classification = classify(cluster_dict_pred, data)
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print(f"{elapsed_time:.2f} minutes")


    data['epitope.pred'] = None

    
    for cluster, tcrs in tqdm(cluster_dict_pred.items()):
        predicted_epitope = cluster_classification[cluster]
        data.loc[data['cdr3'].isin(tcrs), 'epitope.pred'] = predicted_epitope


    ### start ham ###

    couples_ham = find_sequences_within_distance(cdr3, max_mutations)

    start_time = time.time()
    cpm.write_couples_file(couples_ham, "{}/{}".format(out_folder, cdr3_header), "{}_ham".format(out_name))
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print(f"Time taken to write full couples: {elapsed_time:.2f} minutes")

    start_time = time.time()
    cluster_dict_pred = map_clusters(couples_ham, max_neig)
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print(f"{elapsed_time:.2f} minutes")
    # cluster_dict_true = get_true_clusters(data, cdr3_header, epitope_header)

    start_time = time.time()
    cluster_classification = classify(cluster_dict_pred, data)
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print(f"{elapsed_time:.2f} minutes")


    data['epitope.ham'] = None

    
    for cluster, tcrs in tqdm(cluster_dict_pred.items()):
        predicted_epitope = cluster_classification[cluster]
        data.loc[data['cdr3'].isin(tcrs), 'epitope.ham'] = predicted_epitope



    ### end ham ###




    data.to_csv("{}/{}/predicted clusters.csv".format(out_folder, cdr3_header), encoding='utf-8', index=False)

    # data = pd.read_csv("output_files/tests/cdr3/predicted clusters.csv")

    data_chepy = data.dropna(subset=["epitope.pred"])

    # Calculate accuracy
    accuracy = accuracy_score(data_chepy['antigen.epitope'], data_chepy['epitope.pred'])

    # Calculate precision
    precision = precision_score(data_chepy['antigen.epitope'], data_chepy['epitope.pred'], average='weighted', zero_division=0)

    # Calculate recall
    recall = recall_score(data_chepy['antigen.epitope'], data_chepy['epitope.pred'], average='weighted', zero_division=0)

    f1 = f1_score(data_chepy['antigen.epitope'], data_chepy['epitope.pred'], average='weighted', zero_division=0)

    # Print the results
    print(f'chephy Accuracy: {accuracy}')
    print(f'chephy Precision: {precision}')
    print(f'chephy Recall: {recall}')
    print(f'chephy f1_score: {f1}')
    print()


    data_ham = data.dropna(subset=["epitope.ham"])

    # Calculate accuracy
    accuracy = accuracy_score(data_ham['antigen.epitope'], data_ham['epitope.ham'])

    # Calculate precision
    precision = precision_score(data_ham['antigen.epitope'], data_ham['epitope.ham'], average='weighted', zero_division=0)

    # Calculate recall
    recall = recall_score(data_ham['antigen.epitope'], data_ham['epitope.ham'], average='weighted', zero_division=0)

    f1 = f1_score(data_ham['antigen.epitope'], data_ham['epitope.ham'], average='weighted', zero_division=0)

    # Print the results
    print(f'ham Accuracy: {accuracy}')
    print(f'ham Precision: {precision}')
    print(f'ham Recall: {recall}')
    print(f'ham f1_score: {f1}')



if __name__ == "__main__":
    main()