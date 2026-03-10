import numpy as np
import pandas as pd
from collections import defaultdict

# Average Masses of Amino Acids (in Da)
AA_MASS = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841,
    'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406,
    'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
    'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333,
}

def safe_comb(n, k):
    """Compute combination safely (avoiding issues with negative values)"""
    if n < 0 or k < 0 or k > n:
        return 0
    return np.math.comb(n, k)

def b_ion_prob(N, n, i):
    """Probability of b-ion"""
    return safe_comb(N - i - 1, n - 1) / safe_comb(N - 1, n)

def y_ion_prob(N, n, i):
    """Probability of y-ion"""
    return safe_comb(N - i - 1, n - 1) / safe_comb(N - 1, n)

def internal_frag_prob(N, n, length, i, j):
    if n < 2:
        return 0
    elif n == 2:
        return 1 / safe_comb(N - 1, n)
    else:
        return safe_comb(N - 2 - (j - i), n - 2) / safe_comb(N - 1, n)

def aggregate_by_mass(sequence, N, n):
    b_ions = []
    y_ions = []
    internal_mass_counts = defaultdict(float)
    internal_labels_counts = defaultdict(list)

    # b-ions
    for i in range(1, N):
        label = f"b{i}"
        seq = sequence[:i]
        mass = round(sum(AA_MASS[aa] for aa in seq) + 1.00727, 5)
        count = b_ion_prob(N, n, i) * num_simulation
        b_ions.append((label, seq, mass, count))

    # y-ions
    for i in range(1, N):
        label = f"y{i}"
        seq = sequence[-i:]
        mass = round(sum(AA_MASS[aa] for aa in seq) + 18.01560, 5)
        count = y_ion_prob(N, n, i) * num_simulation
        y_ions.append((label, seq, mass, count))

    # Internal fragments
    for length in range(1, N - 1):
        for i in range(1, N - length):
            seq = sequence[i:i + length]
            mass = round(sum(AA_MASS[aa] for aa in seq) + 1.00727, 5)
            probability = internal_frag_prob(N, n, length, i, i + length) * num_simulation
            internal_mass_counts[mass] += probability
            internal_labels_counts[mass].append(f"I{i+1}-{i+length}")

    return b_ions, y_ions, internal_mass_counts, internal_labels_counts

def write_to_excel(b_ions, y_ions, internal_mass_counts, internal_labels_counts):
    """Write the aggregated fragment information to separate sheets in an Excel file."""

    # Convert to DataFrames
    b_ions_df = pd.DataFrame(b_ions, columns=["b-ion Label", "b-ion Sequence", "b-ion Mass", "b-ion Count"])
    y_ions_df = pd.DataFrame(y_ions, columns=["y-ion Label", "y-ion Sequence", "y-ion Mass", "y-ion Count"])

    internal_fragments = []
    for mass, count in internal_mass_counts.items():
        labels = ', '.join(internal_labels_counts[mass])
        internal_fragments.append([labels, mass, count])
    
    internal_fragments_df = pd.DataFrame(internal_fragments, columns=["Internal Fragment Label", "Internal Fragment Mass", "Internal Fragment Count"])

    # Save to Excel with multiple sheets
    with pd.ExcelWriter('Expected_ion_frequencies.xlsx', engine='openpyxl') as writer:
        b_ions_df.to_excel(writer, sheet_name='b-Ions', index=False)
        y_ions_df.to_excel(writer, sheet_name='y-Ions', index=False)
        internal_fragments_df.to_excel(writer, sheet_name='Internal Fragments', index=False)

    print("Results saved to 'Expected_ion_frequencies.xlsx'")

if __name__ == "__main__":
    sequence = input("Please enter the protein sequence: ").upper()

    for aa in sequence:
        if aa not in AA_MASS:
            print(f"Invalid amino acid: {aa}")
            exit()

    N = len(sequence)
    n = int(input("Please enter the number of cleavages: "))
    num_simulation = int(input("Enter the number of simulations you'd like to perform: "))

    b_ions, y_ions, internal_mass_counts, internal_labels_counts = aggregate_by_mass(sequence, N, n)
    write_to_excel(b_ions, y_ions, internal_mass_counts, internal_labels_counts)
