import sympy as sp
import random
import math
from flint import fmpz_mat, fmpq_mat, fmpq, fmpz
from decimal import Decimal, getcontext
import time

class Utils:
    def gram_schmidt(matrix_fmpz):
        matrix = fmpq_mat(matrix_fmpz)
        n, m = matrix.nrows(), matrix.ncols()
        
        result = fmpq_mat(n, m)
        w_array = [[fmpq(0) for _ in range(m)] for _ in range(n)]

        for i in range(n):
            v_i = [matrix[i, j] for j in range(m)]
            w_i = v_i[:]
            
            for j in range(i):
                w_j = w_array[j]
                
                # Calculate dot products
                dot_v_w = sum(v_i[k] * w_j[k] for k in range(m))
                dot_w_w = sum(w_j[k] * w_j[k] for k in range(m))
                
                # Perform subtraction
                factor = dot_v_w / dot_w_w
                w_i = [w_i[k] - factor * w_j[k] for k in range(m)]
            
            w_array[i] = w_i
            
            for j in range(m):
                result[i, j] = w_i[j]

        return result