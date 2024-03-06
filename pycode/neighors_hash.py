import pandas as pd
import itertools
import hashlib
from tqdm import tqdm
import json


def hash_sequences(sequences, hash_method, sequence_map=None):
    sequences_hashed = []
    for seq in sequences:
        hash_value = hash_method(seq)
        if sequence_map is not None:
            sequence_map[hash_value] = seq
        sequences_hashed.append(hash_value)
    return sequences_hashed


df = pd.read_csv('files/atchley.csv')
df.set_index('amino.acid', inplace=True)
atchley_features = df.to_dict(orient='index')
for key, value in atchley_features.items():
    atchley_features[key] = list(value.values())


def hash_atchley_sha_int(sequence, atchley_features=atchley_features):
    hash_input = ''.join([''.join(map(str, atchley_features.get(aa))) for aa in sequence])
    hash_object = hashlib.sha256(hash_input.encode())
    hash_number = int.from_bytes(hash_object.digest(), byteorder='big')
    return hash_number


def generate_mutations(sequence, number):
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    mutations = set()
    for i in range(1, number + 1):
        for positions in itertools.combinations(range(len(sequence)), i):
            for new_aas in itertools.product(aa, repeat=i):
                mutated_sequence = list(sequence)
                for pos, new_aa in zip(positions, new_aas):
                    mutated_sequence[pos] = new_aa
                mutated_sequence = "".join(mutated_sequence)
                if mutated_sequence != sequence:
                    mutations.add(mutated_sequence)
    return mutations


def find_duplicates(lst):
    duplicates = set()
    seen = set()
    for item in lst:
        if item in seen:
            duplicates.add(item)
        else:
            seen.add(item)
    return list(duplicates)


def get_distance(seq1, seq2, distances_df):
    distance = sum([distances_df.loc[aa1, aa2] for aa1, aa2 in zip(seq1, seq2) if aa1 != aa2])
    return distance


sequence_map = {}

sequences_csv = 'files/forchen_F_26L.csv'
sequences_header = "cdr3_amino_acid"

# sequences_csv = 'files/for_chen_B.csv'
# sequences_header = "CDR3.aa"

# sequences_csv = "../random_sequences.csv"
# sequences_header = "sequences"

sequences_df = pd.read_csv(sequences_csv)
sequences_set = set(sequences_df[sequences_header].unique())

sequences_hashed = hash_sequences(sequences_set, hash_atchley_sha_int, sequence_map)


distances_csv = "distance_matrix.csv"
distances_df = pd.read_csv(distances_csv, index_col=0)

max_mutations = 2
couples = {}

for seq in tqdm(sequences_set, desc="Processing {} sequences with {} mutations".format(len(sequences_set), max_mutations)):
    couples[seq] = []
    mutations = generate_mutations(seq, max_mutations)
    mutations_hashed = hash_sequences(mutations, hash_atchley_sha_int)
    combined = sequences_hashed + mutations_hashed
    duplicates = find_duplicates(combined)
    if duplicates:
        for dup in duplicates:
            dup = sequence_map[dup]
            couples[seq].append([dup, get_distance(seq, dup, distances_df)])


couples = {key: sorted(value, key=lambda x: x[1]) for key, value in couples.items()}

output_file = "couples.json"
with open(output_file, "w") as f:
    json.dump(couples, f)
