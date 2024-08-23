import numpy as np
import random
import sympy as sp
import symengine as sym
import math
from flint import fmpz_mat, fmpz, fmpq, fmpq_mat
from decimal import Decimal, getcontext
import time

def random_unimodular(dim, mix):
    T = sp.eye(dim)
    x = sp.zeros(1, dim)
    choices = [-1, 0, 1]
    weights = [1, 5, 1]
    
    for _ in range(mix):
        rows = list(range(dim))
        random.shuffle(rows)
        for i in rows:
            x[0, i] = 1
            for k in range(dim):
                if k != i:
                    x[0, k] = random.choices(choices, weights=weights)[0]
            
            # Perform matrix multiplication more efficiently
            new_row = x * T
            for j in range(dim):
                T[i, j] = new_row[0, j]
            
            # Reset x for the next iteration
            x[0, i] = 0

    return sympy_to_fmpz_mat(T)

def sympy_to_fmpz_mat(basis_sympy):
    return fmpz_mat([[int(item) for item in sublist] for sublist in basis_sympy.tolist()])

def flint_identity_matrix(n):
    identity = []
    for i in range(n):
        row = [0] * n
        row[i] = 1
        identity.append(row)
    return fmpz_mat(identity)



time_now = time.time()
print(random_unimodular(500,2).det())
print(time.time()-time_now)
