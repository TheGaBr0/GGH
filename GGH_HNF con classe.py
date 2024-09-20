from GGH_crypto import GGHHNFCryptosystem, Utils
import sympy as sp
import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpq
import time
import os
import pickle

dimension = 400

GGH_object = GGHHNFCryptosystem(dimension = dimension, GGH_private=False, alpha=0.75, debug=True)
GGH_object.encrypt()

B = GGH_object.public_key[0]

R = GGH_object.private_key

print(fmpz_mat(matrix))

error = GGH_object.error

ciphertext = GGH_object.ciphertext

decrypted_message = GGH_object.decrypt()

print(f"error: {error}")
print(f"error: {decrypted_message}")

n = B.nrows()

# Create the initial matrix with R and a column of zeros
matrix_emb = fmpz_mat([[int(B[i,j]) if j < n else 0 for j in range(n+1)] for i in range(n)])

# Add t as the last row, with 1 as the last element
last_row = [int(ciphertext[0,i]) for i in range(n)] + [1]
matrix_emb = fmpz_mat(matrix_emb.tolist() + [last_row])

time_now = time.time()
found = False
while True:

    matrix_emb, error_message = Utils.BKZ_reduction(matrix_emb, block=20, pruned=True, precision=100, bkzmaxloops=5)

    i_th_row = fmpz_mat([[matrix_emb[0, j] for j in range(n)]])


    if error == i_th_row or error == -1*i_th_row or error_message == None:
        if error == i_th_row or error == -1*i_th_row:
            found = True
            print(f"error: {i_th_row}")
        break
    else:
        print(error_message)

if not found:
    while True:
        matrix_emb, error_message = Utils.BKZ_reduction(matrix_emb, precision=100, block=60, pruned=True, nolll=True, bkzmaxloops=5)

        i_th_row = fmpz_mat([[matrix_emb[0, j] for j in range(n)]])


        if error == i_th_row or error == -1*i_th_row  or error_message == None:
            print(f"error: {i_th_row}")
            break

print(time.time()-time_now)