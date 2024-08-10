from GGH import GGHCryptosystem
from flint import *
import time
import numpy as np
from decimal import Decimal

def row_norm(row):
    return Decimal(sum(Decimal(int(x)**2) for x in row)).sqrt()

def embedding(R, t):
    n = R.nrows()

    # Create the initial matrix with R and a column of zeros
    matrix_emb = fmpz_mat([[int(R[i,j]) if j < n else 0 for j in range(n+1)] for i in range(n)])

    # Add t as the last row, with 1 as the last element
    last_row = [int(t[0,i]) for i in range(n)] + [1]
    matrix_emb = fmpz_mat(matrix_emb.tolist() + [last_row])
    
    matrix_emb = matrix_emb.lll()
    print(matrix_emb)
    min_norm = np.inf
    shortest_row = t

    for i in range(n+1):
        i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n+1)]])
        norm = row_norm(i_th_row)
        if norm < min_norm:
            min_norm = norm
            # We take the first n values of the shortest row
            shortest_row = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])
    
    return t - shortest_row

dimension = 3

R = fmpz_mat([[12, -4, -1], [1,  8, -1], [-4,  1, 14]])
B = fmpz_mat([[145, -73, -23],[-39,  21,  16],[-165,  80,  11]])
m = fmpz_mat([[-48, 29, -76]])
e = fmpz_mat([[3, 3, 3]])

start_time = time.time()
GGH_object = GGHCryptosystem(dimension = dimension, R=R, B=B, m=m, e=e)
GGH_object.encrypt()

R = GGH_object.private_basis
B = GGH_object.public_basis
U = GGH_object.unimodular

print(R)
print(f"R: {GGH_object.get_hadamard_ratio(R)}")
print(B)
print(f"B: {GGH_object.get_hadamard_ratio(B)}")
print(U)


message = GGH_object.message

ciphertext = GGH_object.ciphertext

print(f"C: {ciphertext}")

decrypted_message = GGH_object.decrypt()

print(f"time taken: {time.time() - start_time}")

print(f"sigma: {GGH_object.public_key[1]}")
print(f"error: {GGH_object.error}")
print(f"message: {message}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"decrypted_message: {decrypted_message}")
print(decrypted_message == message)
print(embedding(B, ciphertext))
print(embedding(B, ciphertext) * B.inv())
        
a = ciphertext * B.inv()
print(a)
cols = a.ncols()
rows = a.nrows()

for i in range(rows):
    for j in range(cols):
        a[i,j] = round(a[i,j])

print(a)