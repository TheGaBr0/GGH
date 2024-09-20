from GGH_crypto import Utils, GGHCryptosystem
from flint import fmpz_mat
from decimal import Decimal
import time
import os
import pickle
import lzma

dimension = 400

GGH = GGHCryptosystem(dimension = dimension, debug=True)
GGH.encrypt()

def write_and_get_size(matrix, filename, use_pickle=False):
    if use_pickle:
        with lzma.open(filename, 'wb') as f:
            pickle.dump(matrix.tolist(), f)
    else:
        Utils.write_matrix_to_file(matrix, filename)
    return os.path.getsize(filename) / (1024 * 1024)

# Calculate sizes for text files
R_txt = write_and_get_size(GGH.private_key, "priv.txt")
B_txt = write_and_get_size(GGH.public_basis, "pub.txt")
C_txt = write_and_get_size(GGH.ciphertext, "cip.txt")

# Calculate sizes for pickle files
R_pkl = write_and_get_size(GGH.private_key, "priv.pkl.xz", use_pickle=True)
B_pkl = write_and_get_size(GGH.public_basis, "pub.pkl.xz", use_pickle=True)
C_pkl = write_and_get_size(GGH.ciphertext, "cip.pkl.xz", use_pickle=True)

# Print results
print("Text file sizes (MB):")
print(f"Private key: {R_txt:.2f}")
print(f"Public basis: {B_txt:.2f}")
print(f"Ciphertext: {C_txt:.2f}")

print("\nPickle file sizes (MB):")
print(f"Private key: {R_pkl:.2f}")
print(f"Public basis: {B_pkl:.2f}")
print(f"Ciphertext: {C_pkl:.2f}")


