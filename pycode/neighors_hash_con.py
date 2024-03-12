import pandas as pd
import itertools
import hashlib
from tqdm import tqdm
import json

df = pd.read_csv('files/atchley.csv')
df.set_index('amino.acid', inplace=True)
atchley_features = df.to_dict(orient='index')


def hash_atchley_sha_int(sequence):
    hash_input = tuple(atchley_features.get(aa) for aa in sequence)
    hash_object = hashlib.sha256(str(hash_input).encode())
    hash_number = int.from_bytes(hash_object.digest(), byteorder='big')
    return hash_number


def hash_sequences(mutations_lst, hash_method):
    new_lst = []
    for row in tqdm(mutations_lst, desc="hashing sequences"):
        hash_value = hash_method(row[0])
        new_lst.append([row[0], hash_value, row[2]])
    return new_lst


def generate_mutations(sequence, number, mutations_lst):
    mutations_lst.append([sequence, 0, ""])
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    for i in range(1, number + 1):
        for positions in itertools.combinations(range(len(sequence)), i):
            for new_aas in itertools.product(aa, repeat=i):
                mutated_sequence = list(sequence)
                for pos, new_aa in zip(positions, new_aas):
                    mutated_sequence[pos] = new_aa
                mutated_sequence = "".join(mutated_sequence)
                if mutated_sequence != sequence:
                    mutations_lst.append([mutated_sequence, 0, sequence])
                    # mutations_df = mutations_df.append({'mutation': mutated_sequence, 'original': sequence}, ignore_index=True)

    return mutations_lst


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

mutations_df = pd.DataFrame(columns=['mutation','mutation_hash' 'original'])

sequences_df = pd.read_csv(sequences_csv)
sequences_set = set(sequences_df[sequences_header].unique())
# sequences_hashed = hash_sequences(sequences_set, hash_atchley_sha_int,sequence_map)


max_mutations = 1

couples = {}
# pd.DataFrame(columns=['mutation','mutation_hash' 'original'])
mutations_lst = []


for seq in tqdm(sequences_set, desc="generating mutations"):
    couples[seq] = []
    mutations_lst = generate_mutations(seq, max_mutations, mutations_lst)


mutations_lst = hash_sequences(mutations_lst, hash_atchley_sha_int)

mutations_df = pd.DataFrame(mutations_lst, columns=['mut', 'mut_hash', 'org'])
mutations_sorted = mutations_df.sort_values(by='mut_hash')


distances_csv = "distance_matrix.csv"
distances_df = pd.read_csv(distances_csv, index_col=0)

couples = {}

# for seq in tqdm(sequences_set, "finding close sequences"):

# couples[seq] = []
# filtered_df = mutations_sorted[(mutations_sorted['org'] == seq)\
#                                 | (mutations_sorted['org'] == "")]

duplicate_values = mutations_sorted[mutations_sorted.duplicated(['mut_hash'])]

for index, row in duplicate_values.iterrows():
    seq = row["org"]
    dup = row["mut"]
    if seq != "":
        if not seq in couples:
            couples[seq] = []
        # print(row["mut"], row["org"])   
        couples[seq].append([dup, get_distance(seq, dup, distances_df)])


couples = {key: sorted(value, key=lambda x: x[1]) for key, value in couples.items()}

output_file = "couples.json"
with open(output_file, "w") as f:
    json.dump(couples, f)
