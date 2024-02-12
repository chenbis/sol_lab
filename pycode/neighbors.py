import json
import pandas as pd

# sorts amino acids couples according to the distance df

def get_sorted_substitutions(distances_df):
    substitutions = []
    amino_acids = distances_df.index

    for source_acid in amino_acids:
        for target_acid in amino_acids:
            if source_acid != target_acid:
                distance = distances_df.loc[source_acid, target_acid]
                substitutions.append((source_acid, target_acid, distance))

    return sorted(substitutions, key=lambda x: x[2])


def find_couples(sequences_set, sorted_amino, max_neighbors=1):
    couples = {}
    for seq in sequences_set:
        letter_indexes = {}
        couples[seq] = set()
        flag = 0
        for sub in sorted_amino:
            
            if flag < max_neighbors:
                if not sub[0] in letter_indexes:
                    letter_indexes[sub[0]] = [index for index, char in enumerate(seq) if char == sub[0]]
                for occ in letter_indexes[sub[0]]:
                    seq_to_search = seq[:occ] + sub[1] + seq[occ + 1:]
                    # if seq_to_search in couples:
                    #     continue
                    if seq_to_search in sequences_set:
                        
                        # print(seq, seq_to_search)
                        couples[seq].add(seq_to_search)
                        flag += 1
                        if flag < max_neighbors:
                            break

                
            else:
                continue

    return couples


# read distance matrix
distances_csv = "distance_matrix.csv"

sequences_csv = 'files/forchen_F_26L.csv'
sequences_header = "cdr3_amino_acid"

# sequences_csv = 'files/for_chen_B.csv'
# sequences_header = "CDR3.aa"

# sequences_csv = "../random_sequences.csv"
# sequences_header = "sequences"


distances_df = pd.read_csv(distances_csv, index_col=0)

sorted_amino = get_sorted_substitutions(distances_df)

sequences_df = pd.read_csv(sequences_csv)

sequences = sequences_df[sequences_header].tolist()
sequences_set = set(sequences)

# sequences_set.add("CASSLALAGGTDTQYI")
# sequences_set.add("CASSLALAGGTDTQYC")

sorted_amino = get_sorted_substitutions(distances_df)
couples = find_couples(sequences_set, sorted_amino, max_neighbors=3)

couples_list = {key: list(value) for key, value in couples.items()}

# Serialize data into file without escape characters:
with open("couples.json", 'w') as outfile:
    json.dump(couples_list, outfile, indent=4)