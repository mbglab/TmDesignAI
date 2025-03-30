# Copyright (C) 2022-2025, Bin-Guang Ma (mbg@mail.hzau.edu.cn). All rights reserved.
# This script is for the E. coli Terminator Strength Prediction Model.
# Author: Jie Li, 2022-12; Revised by Lin-Feng Wu and Bin-Guang Ma, 2025-03. 

import numpy as np
import itertools
from scipy import stats

# Load file and shuffle data
file_path = 'shujudaquan_new.txt'
with open(file_path, 'r') as f:
    data = f.readlines()[1:]  # Skip the header line

np.random.seed(1314)
np.random.shuffle(data)


# Convert Uracil (U) to Thymine (T)
def convert_U_to_T(seq):
    """Replace 'U' with 'T' in the given sequence."""
    return seq.replace('U', 'T')


# Convert lowercase letters to uppercase
def to_uppercase(seq):
    """Convert all lowercase letters in a sequence to uppercase."""
    return seq.upper()


# Calculate overlap length between two sequences
def calculate_overlap(seq1, seq2):
    """
    Returns the overlap length between two sequences.
    If two sequences are identical, returns -1.
    """
    min_length = min(len(seq1), len(seq2))
    overlap_length = 0

    for i in range(min_length):
        if seq1[-i-1:] == seq2[:i+1]:
            overlap_length = i + 1

    return -1 if seq1 == seq2 else overlap_length


# Merge two sequences based on overlap
types = []  # Stores the type of merge operation
def merge_sequences(seq1, seq2):
    """
    Merge two sequences based on their overlap.
    - If one sequence contains the other, return the longer sequence.
    - If there is a sufficient overlap (>20 bp), merge them accordingly.
    - Otherwise, return 2 as an indicator of insufficient overlap.
    """
    overlap1 = calculate_overlap(seq1, seq2)
    overlap2 = calculate_overlap(seq2, seq1)

    if seq1 in seq2:
        types.append(0)
        return seq2
    if seq2 in seq1:
        types.append(0)
        return seq1
    if overlap1 > 20:
        types.append(-1)
        return seq1 + seq2[overlap1:]
    if overlap2 > 20:
        types.append(1)
        return seq2 + seq1[overlap2:]

    types.append(2)
    return 2

# Extract sequences
sequences = []
for line in data:
    tmp = line.strip().split('\t')
    seq1 = to_uppercase(tmp[0])  # Convert to uppercase
    seq2 = convert_U_to_T(tmp[2]) + convert_U_to_T(tmp[3]) + convert_U_to_T(tmp[5])
    #sequences.append(seq_union(seq1,seq2))
    sequences.append(seq2)

# data = [data[i] for i in range(len(types)) if types[i] != 2]
# sequences = [seq for seq in sequences if seq != 2]

# Count occurrences of a substring in a given sequence
def count_substring_occurrences(seq, substring):
    """
    Count the occurrences of a specific substring in a given sequence.
    """
    count = 0
    seq_len = len(seq)
    sub_len = len(substring)

    for i in range(seq_len - sub_len + 1):
        if seq[i:i + sub_len] == substring:
            count += 1
    return count


# Composition of k-spaced Nucleic Acid Pairs (CKSNAP)
def encode_CKSNAP(seq, k_spaced=1, remove_zeros=True):
    """
    Compute the CKSNAP (Composition of k-spaced Nucleic Acid Pairs) feature representation.

    Parameters:
    - seq (str): The input DNA/RNA sequence.
    - k_spaced (int): The spacing between nucleotides.
    - remove_zeros (bool): If True, remove the first 16 elements from the feature vector.

    Returns:
    - List of CKSNAP features.
    """
    CKSNAP = []
    di_nucleotide_dict = {}

    for i in range(k_spaced + 1):
        # Initialize all possible dinucleotide pairs (AA, AC, ..., TT)
        for pair in itertools.product("ACGT", repeat=2):
            di_nucleotide_dict[''.join(pair)] = 0

        # Count occurrences of spaced dinucleotide pairs
        for k in range(len(seq) - i - 1):
            di_nucleotide_dict[''.join((seq[k], seq[k + i + 1]))] += 1

        # Normalize frequencies
        for key in di_nucleotide_dict.keys():
            di_nucleotide_dict[key] /= (len(seq) - i - 1)

        # Append to feature list
        CKSNAP.extend([di_nucleotide_dict[''.join(pair)] for pair in itertools.product("ACGT", repeat=2)])

    # Optionally remove the first 16 elements
    return CKSNAP[16:] if remove_zeros else CKSNAP


# K-mer Encoding
def encode_Kmer(seq, K=2, include_lower_orders=True):
    """
    Compute K-mer frequency representation for the given sequence.

    Parameters:
    - seq (str): The input DNA/RNA sequence.
    - K (int): The length of K-mers.
    - include_lower_orders (bool): If True, also compute features for all k-mers from 1 to K.

    Returns:
    - List of K-mer frequencies.
    """
    assert K >= 1, 'K should be at least 1!'
    kmer_frequencies = []

    if not include_lower_orders:
        for kmer in itertools.product("ACGT", repeat=K):
            kmer_frequencies.append(count_substring_occurrences(seq, ''.join(kmer)) / (len(seq) - K + 1))
    else:
        for j in range(1, K + 1):
            for kmer in itertools.product("ACGT", repeat=j):
                kmer_frequencies.append(count_substring_occurrences(seq, ''.join(kmer)) / (len(seq) - j + 1))

    return kmer_frequencies


# Read Pseudo Nucleic Acid Composition (PseNAC) features from file
def read_PseNAC(file_path):
    """
    Read Pseudo Nucleic Acid Composition (PseNAC) features from a file.

    Parameters:
    - file_path (str): Path to the PseNAC feature file.

    Returns:
    - List of feature vectors.
    """
    features = []
    with open(file_path, 'r') as file:
        for line in file:
            features.append([float(value.split(':')[1]) for value in line[1:].strip().split(' ') if value])

    return features


# Extract sequence features using 3-mer encoding
sequence_features = [encode_Kmer(seq, K=3, include_lower_orders=False) for seq in sequences]
print(len(sequence_features[0]))  # Print the length of one feature vector


# === Extract Thermodynamic Features ===

# Extract features from 2013 Nucleic Acids Research (NAR) study
NAR2013_features = []
ULHAB_features = []

for line in data:
    tmp = line.strip().split('\t')
    ULHAB_features.append([float(tmp[6]), float(tmp[7]), float(tmp[8]), float(tmp[10])])


# === Extract Minimum Free Energy (MFE) Features ===
mfe_file_path = 'seq606_RNAfold.txt'

with open(mfe_file_path, 'r') as mfe_file:
    data_MFE = mfe_file.readlines()

MFE_features = [
    float(line.strip().split(' (')[1][:-1]) for line in data_MFE if '-' in line
]


# === Thermodynamic Calculation Functions ===
def calculate_stem_length(hairpin, loop):
    """Compute the length of the Terminator stem region."""
    return (len(hairpin) - len(loop)) // 2


def calculate_GH_L(GH, length):
    """Compute the normalized GH feature by sequence length."""
    return GH / length


def calculate_GH_L2(GH, stack_length):
    """Compute the normalized GH feature by stack length."""
    return GH / stack_length


# === Merge Thermodynamic Features ===
features_relixue = []
for i, line in enumerate(data):
    tmp = line.strip().split('\t')
    features_relixue.append(
        ULHAB_features[i] +
        [MFE_features[i]] +
        [calculate_GH_L(float(tmp[8]), len(tmp[3]))] +
        [calculate_GH_L2(float(tmp[8]), calculate_stem_length(tmp[3], tmp[4]))]
    )


# === Tail Region Feature Calculation ===
def compute_tail_score(seq):

    u = [1] + [0] * len(seq)
    for i in range(len(seq)):
        u[i + 1] = u[i] * (0.9 if seq[i] in 'TU' else 0.6)
    return sum(u[1:])


def compute_closing_stack(seq):

    if seq[0] == 'G':
        return 2
    elif seq[0] == 'C':
        return 1
    elif seq[0] in 'UT':
        return 0
    else:
        return -1


# Compute NAR2013-based features
NAR2013_features = [
    [compute_tail_score(tmp[5]), compute_tail_score(tmp[5][:5]), compute_tail_score(tmp[5][5:8]), compute_closing_stack(tmp[3])]
    for tmp in (line.split('\t') for line in data)
]


# === Extract Features Based on Secondary Structure ===
new2022_features = [
    [calculate_stem_length(tmp[3], tmp[4])]
    for tmp in (line.split('\t') for line in data)
]


# === Sequence Feature Extraction Based on Secondary Structure Regions ===
def extract_region_features(seq, mode):
    """
    Extract sequence features based on different secondary structure regions.
    - mode 1: CKSNAP(k=2) + Kmer(k=2)
    - mode 2: Kmer(k=2) only
    - mode 3: CKSNAP(k=3) + Kmer(k=3)
    """
    if mode == 1:
        return encode_CKSNAP(seq, k_spaced=2) + encode_Kmer(seq, K=2)
    elif mode == 2:
        return encode_Kmer(seq, K=2)
    elif mode == 3:
        return encode_CKSNAP(seq, k_spaced=3) + encode_Kmer(seq, K=3)
    else:
        raise ValueError("Invalid mode! Choose from {1, 2, 3}.")


# Extract features for different RNA regions
fenquyu_features = [
    extract_region_features(convert_U_to_T(to_uppercase(tmp[2])), mode=1) +
    extract_region_features(convert_U_to_T(to_uppercase(tmp[3])), mode=1) +
    extract_region_features(convert_U_to_T(to_uppercase(tmp[4])), mode=2) +
    extract_region_features(convert_U_to_T(to_uppercase(tmp[5])), mode=1)
    for tmp in (line.split('\t') for line in data)
]

# Alternative mode 3 (if needed)
# fenquyu_features = [
#     extract_region_features(convert_U_to_T(to_uppercase(tmp[2])), mode=3) +
#     extract_region_features(convert_U_to_T(to_uppercase(tmp[3])), mode=3) +
#     extract_region_features(convert_U_to_T(to_uppercase(tmp[4])), mode=2) +
#     extract_region_features(convert_U_to_T(to_uppercase(tmp[5])), mode=3)
#     for tmp in (line.split('\t') for line in data)
# ]


# === Merge All Features ===
features = [
    sequence_features[i] +
    features_relixue[i] +
    NAR2013_features[i] +
    new2022_features[i] +
    fenquyu_features[i]
    for i in range(len(data))
]

import numpy as np
import math
import pandas as pd
from sklearn import metrics
from sklearn.feature_selection import f_regression, SelectKBest
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn import svm

# === Extract Labels ===
values = []
pre_2013NM = []

for line in data:
    cols = line.strip().split('\t')
    values.append(math.log(float(cols[1]), 2))  # Log transformation for stability
    pre_2013NM.append(math.log(float(cols[-1]), 2))  # Log transformation for predictions

# Convert to NumPy arrays
features = np.array(features)
values = np.array(values)

# === Split Dataset into Training and Testing Sets ===
split_index = len(features) // 5  # Use 80% training, 20% testing
train_data, test_data = features[split_index:], features[:split_index]
train_values, test_values = values[split_index:], values[:split_index]

# === Feature Normalization ===
scaler = StandardScaler()
train_data = scaler.fit_transform(train_data)
test_data = scaler.transform(test_data)

# === Feature Selection using f_regression ===
print(f"Number of features: {train_data.shape[1]}")

selector = SelectKBest(score_func=f_regression, k=train_data.shape[1])
results = selector.fit(train_data, train_values)

# Handle NaN values (set them to -1)
PX = np.nan_to_num(results.scores_, nan=-1)

# Rank features by importance
feature_ranking = np.argsort(PX)[::-1]

# === Train SVR  ===
cv_folds = 5  # Number of folds for cross-validation
score_methods = ['r2', 'neg_mean_absolute_error', 'neg_mean_squared_error']
best_scores = []

for score_method in score_methods:
    for num_features in range(100, train_data.shape[1] + 1):
        selected_train_data = train_data[:, feature_ranking[:num_features]]

        parameters = {
            "kernel": ["rbf"],
            "C": 2.0 ** np.arange(-2, 8),
            "gamma": 2.0 ** np.arange(-13, 5)
        }

        grid = GridSearchCV(svm.SVR(), parameters, cv=cv_folds, scoring=score_method, verbose=1)
        grid.fit(selected_train_data, train_values)

        print(
            f"{num_features} features - SVC Best {score_method}: {grid.best_score_}, Best Params: {grid.best_params_}")

        # Store the best results
        best_scores.append((num_features, score_method, grid.best_score_, grid.best_params_))
