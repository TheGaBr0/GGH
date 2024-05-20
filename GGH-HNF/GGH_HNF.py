import numpy as np
import random
import math
from flint import fmpz_mat, fmpq_mat, fmpz
import sympy as sp
from decimal import Decimal, getcontext

class GGHCryptosystem:
    def __init__(self, dimension, use_random=True, integer_sigma=True):
        self.dimension = dimension
        self.integer_sigma = integer_sigma
        
        self.ciphertext = None
        self.error = None

        self.good_basis = None
        self.bad_basis = None

        self.private_key = None
        self.public_key = None

        if use_random:
            self.generate_random_error()
            self.generate_keys()

    def numpy_to_fmpz_mat(self, numpy_matrix):
        return fmpz_mat([[int(item) for item in sublist] for sublist in numpy_matrix.tolist()])

    def column_norm(self, col):
        getcontext().prec = 50
        return Decimal(sum(Decimal(int(x))**2 for x in col)).sqrt()

    def min_norm_column(self, matrix):
        norms = []
        for j in range(matrix.ncols()):
            column = [matrix[i, j] for i in range(matrix.nrows())]
            norms.append(self.column_norm(column))
        min_norm = min(norms)
        return norms[norms.index(min_norm)]

    def generate_rho(self, R): #GGH originale pagina 9, Theorem 1

        min_norm = self.min_norm_column(R)

        rho = Decimal(0.5) * min_norm

        return rho
    
    def generate_random_error(self):
        n = self.dimension
        random_elements = [[random.randint(-n, n) for _ in range(n)]] #pagina 3 "GGH - altre spiegazioni"
        self.error = fmpz_mat(random_elements).transpose()
    
    def reduce_mod_B(self):
        
        r = fmpq_mat(self.error)
        
        H = self.bad_basis
        
        x = fmpz_mat(r.nrows(), r.ncols())
        
        # Iterate over each coordinate starting from the last one
        for i in reversed(range(r.nrows())):
            sum_j = sum(H[i, j] * x[j, 0] for j in range(i + 1, self.dimension))
            # Compute x[i] using the formula
            x[i, 0] = math.floor((r[i, 0] - sum_j) / H[i, i])

        return x

    def generate_keys(self):
        n = self.dimension
        l = n
        k = fmpz(l * math.ceil(math.sqrt(self.dimension) + 1))
        print("Generating private key...")
        while True:
            try:
                R = self.numpy_to_fmpz_mat(np.random.randint(-l, l-1, (self.dimension, self.dimension)))
                I = self.numpy_to_fmpz_mat(np.eye(self.dimension))
                KI = k * I
                R += KI
                R_inv = R.inv()
            except Exception as e:
                print(e)
                continue
            else:
                break
        self.good_basis = R
        rho = self.generate_rho(R)
        print("Generating public key...")

        self.bad_basis = H = R.hnf().transpose()

        self.public_key = H

        self.private_key = (R_inv, R)

    def encrypt(self):
        x = self.reduce_mod_B()
        
        H = self.bad_basis
        r = self.error

        self.ciphertext =  r - H * x

    def decrypt(self):
        R_inv, R = self.private_key
        
        c = self.ciphertext
        x = R_inv.transpose() * c

        
        rounded_x = fmpz_mat(x.nrows(), x.ncols())

        for i in range(x.nrows()):
            for j in range(x.ncols()):
                rounded_x[i, j] = round(x[i,j])
                
        result = R.transpose() * rounded_x
        
        return c - result
    
    def get_hadamard_ratio(self, basis = None, type='public'):
        
        if not basis: 
            if basis not in ['public', 'private']:
                print("Valid inputs are: 'public', 'private'")
                return None
            matrix = self.bad_basis if basis == 'public' else self.good_basis
        else:
            matrix = basis
        norms = []
        
        for j in range(matrix.ncols()):
            column = [matrix[i, j] for i in range(matrix.nrows())]
            norms.append(self.column_norm(column))
        
        denominator = math.prod(norms)
        numerator = abs(Decimal(matrix.det().str()))
        result = (numerator / denominator) ** Decimal(1 / self.dimension)
        return f"{result:.16f}"
