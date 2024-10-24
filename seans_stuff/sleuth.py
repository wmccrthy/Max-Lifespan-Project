# %%
# import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

# %%
file_path = "value_counts.txt"

with open(file_path, 'r') as file:
    lines = file.readlines()
    print(lines)


# %%
file_path = 'sequence_lengths.txt'

# Initialize an empty list to store lengths
sequence_lengths = []

# Read the file in chunks to avoid memory issues
with open(file_path, 'r') as file:
    for line in file:
        length = int(line.strip())  # Convert each line to an integer (assuming they are integer lengths)
        logl = np.log10(length)
        if  (logl > 2) & (logl < 4.5):
            sequence_lengths.append(length)

# Convert list to numpy array for faster processing
sequence_lengths = np.array(sequence_lengths)
min_len, max_len = np.min(sequence_lengths), np.max(sequence_lengths)
print(min_len, max_len)

# Generate the histogram
plt.figure(figsize=(10, 6))
plt.hist(sequence_lengths, bins=50, edgecolor='black')  # Adjust bins for better resolution
plt.title('Histogram of Sequence Lengths')
plt.xlabel('Sequence Length')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

# %%
pth =  "/data/rbg/users/wmccrthy/rbgquanta1/wmccrthy/cumulativeTOGAset_trimmed.csv"
df = pd.read_csv(pth, chunksize=100)
first = next(df)

# %%
# Open the HDF5 file in read mode
file_path = '/data/rbg/users/seanmurphy/maxls/maggies_stuff/prelim_training/data/tokenized_enformer_829124_training.h5'
with h5py.File(file_path, 'r') as f:
    dataset = f['dataset']
    
    # Read the first 3 elements
    first_three_elements = dataset[:3]
    print(first_three_elements)
# %%
