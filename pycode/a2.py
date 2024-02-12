import random, csv

# Amino acid symbols
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

# Function to generate a random amino acid sequence of given length
def generate_random_sequence(length):
    return ''.join(random.choice(amino_acids) for _ in range(length))

# Function to create a sequence with a change of one random amino acid
def mutate_sequence(sequence):
    num_mutations = random.randint(1, int(len(sequence)/2))
    for i in range(num_mutations):
        index = random.randint(0, len(sequence) - 1)
        new_amino_acid = random.choice(amino_acids.replace(sequence[index], ''))
        sequence = sequence[:index] + new_amino_acid + sequence[index + 1:]
    return sequence

# Write output to CSV
def write_to_csv(csv_filename, sequences):
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["sequences"])
        for sequence in sequences:
            writer.writerow([sequence])


sequences = []

for i in range(15):


    # Number of sequences to generate
    num_sequences = random.randint(16,1000)

    # Number of mutations to perform on each sequence
    num_mutations = random.randint(15,1000)


    # Generate initial random sequence
    initial_sequence = generate_random_sequence(random.randint(15,20))
    sequences.append(initial_sequence)

    for i in range(random.randint(1,4)):
        # Generate and output mutated sequences
        for i in range(num_sequences - 1):  # We already have one initial sequence
            mutated_sequence = initial_sequence
            for _ in range(num_mutations):
                mutated_sequence = mutate_sequence(mutated_sequence)
                sequences.append(mutated_sequence)






# Write output to CSV
write_to_csv('random_sequences.csv', sequences)

print("Output saved")
