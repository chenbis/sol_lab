import pandas as pd
from collections import defaultdict
import json, csv, time
from tqdm import tqdm
from Levenshtein import distance
from multiprocessing import Pool

MUTATIONS_DIR = "mutations"
max_mutations = 3

def find_strings_within_distance(strings_set, new_string, max_distance):
    return {string for string in strings_set if string != new_string and distance(string, new_string) <= max_distance}

def truncate_sequences(args):
    csv_file, column_name, right, left = args
    sequences_set = set()
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            sequence = row[column_name].strip()
            start_index = max(0, len(sequence)//2 - left)
            end_index = min(len(sequence), len(sequence)//2 + right)
            sequences_set.add(sequence[start_index:end_index])
    return sequences_set

def process_sequence(args):
    seq, sequences_set = args
    return [(var, seq) for var in find_strings_within_distance(sequences_set, seq, max_mutations)]



def main():

    # 1000 sequences
    sequences_file = 'files/forchen_F_26L.csv'
    sequences_header = "cdr3_amino_acid"

    # # 40K sequences
    # sequences_file = 'files/for_chen_B.csv'
    # sequences_header = "CDR3.aa"

    # # 1M sequences
    # sequences_file = "random_sequences.csv"
    # sequences_header = "sequences"

    sequences_file = 'files/articles/mels.csv'
    sequences_header = "mel7.cdr3a"

    sequences_set = truncate_sequences((sequences_file, sequences_header, 4, 4))

    print("creating mutations")

    with Pool() as pool:
        results = list(tqdm(pool.imap(process_sequence, [(seq, sequences_set) for seq in sequences_set]), total=len(sequences_set)))

    neighbors = {item for sublist in results for item in sublist}


    distances_csv = "distance_matrix.csv"
    distances_df = pd.read_csv(distances_csv, index_col=0)
    distances_dict = {(aa1, aa2): distances_df.loc[aa1, aa2] for aa1 in distances_df.index for aa2 in distances_df.columns}
    
    couples = defaultdict(list)

    for mut, seq in neighbors:
        couples[seq].append([mut, sum(distances_dict[(aa1, aa2)] for aa1, aa2 in zip(seq, mut))])

    couples = {key: sorted(value, key=lambda x: x[1]) for key, value in couples.items()}

    output_file = "output_files/compare/couples_by_che_{}.json".format("test")
    with open(output_file, "w") as f:
        json.dump(couples, f)

    
    cutoff = 1  # Maximum number of changes allowed
    couples = defaultdict(list)
    for seq in sequences_set:
        close_matches = [word for word in sequences_set if word != seq and distance(seq, word) <= cutoff]
        for mut in close_matches:
            couples[seq].append([mut,distance(seq, mut)]) 

    output_file = "output_files/compare/couples_by_lev_{}.json".format("test")
    with open(output_file, "w") as f:
        json.dump(couples, f)


if __name__ == "__main__":
    main()