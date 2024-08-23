from GGH_crypto import GGHCryptosystem
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

dimension = 100


start_time = time.time()
GGH_object = GGHCryptosystem(dimension = dimension)
GGH_object.encrypt()
print(GGH_object.unimodular.det())


message = GGH_object.message
decrypted_message = GGH_object.decrypt()

print(f"message: {message}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"decrypted_message: {decrypted_message}")
print(decrypted_message == message)