from GGH_crypto import GGHHNFCryptosystem, Utils
import sympy as sp
import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpq

dimension = 300

R = fmpz_mat([[12, -4, -1], [1,  8, -1], [-4,  1, 14]])
e = fmpz_mat([[1,1,2]])

GGH_object = GGHHNFCryptosystem(dimension = dimension, GGH_private=False, debug=True)
GGH_object.encrypt()

B = GGH_object.public_key[0]

R = GGH_object.private_key
    
error = GGH_object.error

ciphertext = GGH_object.ciphertext

decrypted_message = GGH_object.decrypt()

print(f"error: {error}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"error: {decrypted_message}")

n = B.nrows()

exit()

# Create the initial matrix with R and a column of zeros
matrix_emb = fmpz_mat([[int(B[i,j]) if j < n else 0 for j in range(n+1)] for i in range(n)])

# Add t as the last row, with 1 as the last element
last_row = [int(ciphertext[0,i]) for i in range(n)] + [1]
matrix_emb = fmpz_mat(matrix_emb.tolist() + [last_row])

matrix_20 = Utils.BKZ_reduction(matrix_emb)

min_norm = float('inf')

for i in range(n + 1):
    i_th_row = fmpz_mat([[matrix_20[i, j] for j in range(n + 1)]])
    norm = Utils.vector_l2_norm(i_th_row)
    if norm < min_norm:
        min_norm = norm
        # Prendi i primi n valori della riga più corta
        shortest_vector = fmpz_mat([[matrix_20[i, j] for j in range(n)]])

print(f"error: {shortest_vector}")

matrix_60 = Utils.BKZ_reduction(matrix_20, block=60, pruned=True, nolll=True)

min_norm = float('inf')
for i in range(n + 1):
    i_th_row = fmpz_mat([[matrix_20[i, j] for j in range(n + 1)]])
    norm = Utils.vector_l2_norm(i_th_row)
    if norm < min_norm:
        min_norm = norm
        # Prendi i primi n valori della riga più corta
        shortest_vector = fmpz_mat([[matrix_20[i, j] for j in range(n)]])

print(f"error: {shortest_vector}")