import numpy as np
import random
import math
from flint import fmpz_mat, fmpq_mat, fmpq, fmpz
from decimal import Decimal, getcontext
import time

class GGHHNFCryptosystem:
    def __init__(self, dimension, R=None, B=None, x=None, e=None, rho_check=True, random_private=True, debug=True):
        self.dimension = dimension
        
        self.x = x
        self.ciphertext = None
        self.error = e

        self.good_basis = R
        self.bad_basis = B
        self.unimodular = None
        
        self.R_rho = None
        self.rho_check = rho_check
        self.random_private = random_private
        self.debug = True

        self.private_key = None
        self.public_key = None

        if R is not None:
            # Check if R is a 2D array or matrix
            if R.nrows() != dimension:
                raise ValueError(f"R must be a {dimension}x{dimension} matrix, but got {R.nrows()}")
            
            if B is None:
                self.generate_keys_from_R()
            else:
                self.bad_basis = B
                self.generate_keys_with_B()
        else:
            self.generate_keys()
            
        if e is None:
            self.generate_random_error()

    def generate_keys_from_R(self):
        R = self.good_basis
        R_inv = R.inv()
        print("Generating public key...")
        self.bad_basis = H = R.hnf()
        print("Computing results...")
        self.public_key = H
        self.private_key = (R_inv, R)

    def generate_keys_with_B(self):
        R = self.good_basis
        R_inv = R.inv()
        print("Computing results...")
        self.public_key = self.bad_basis
        self.private_key = (R_inv, R)

    def numpy_to_fmpz_mat(self, numpy_matrix):
        return fmpz_mat([[int(item) for item in sublist] for sublist in numpy_matrix.tolist()])

    def vector_norm(self, row):
        
        if isinstance(row, fmpz_mat):
            row = fmpq_mat(row)
        getcontext().prec = 50
        return Decimal(sum((Decimal(int(x.numer())) / Decimal(int(x.denom()))) ** 2 for x in row)).sqrt()

    def min_norm_row(self, matrix):
        norms = []
        for j in range(matrix.nrows()):
            row = [matrix[j, i] for i in range(matrix.ncols())]
            norms.append(self.vector_norm(row))
        min_norm = min(norms)
        return norms[norms.index(min_norm)]

    def calculate_rho(self, basis):
        basis_orthogonalized = gram_schmidt(basis)
        min_norm = self.min_norm_row(basis_orthogonalized)
        rho = Decimal(0.5) * min_norm
        return float(f'{rho:f}')
    
    def generate_random_error(self):
        n = self.dimension
        if self.rho_check:
            print("Calculating rho...")
            R_rho = self.calculate_rho(self.good_basis)
            rand_val = n
            
            random_elements = [[random.randint(-rand_val, rand_val) for _ in range(n)]]
            error = fmpz_mat(random_elements)
            
            while(not self.vector_norm(error) < R_rho):
                rand_val -= 1
                random_elements = [[random.randint(-rand_val, rand_val) for _ in range(n)]]
                error = fmpz_mat(random_elements)
            
            self.error = error
        else:
            random_elements = [[random.randint(-28, 28) for _ in range(n)]]
            self.error = fmpz_mat(random_elements)
    
    def reduce_mod_B(self):
        r = fmpq_mat(self.error)
        H = self.bad_basis
        x = fmpz_mat(r.nrows(), r.ncols())
        
        # Iterate over each coordinate starting from the last one
        for i in reversed(range(r.ncols())):
            sum_j = sum(H[j, i] * x[0, j] for j in range(i + 1, self.dimension))
            # Compute x[i] using the formula
            x[0, i] = math.floor((r[0, i] - sum_j) / H[i, i])

        return x

    def generate_keys(self):
        n = self.dimension
        if not self.random_private:
            time_start = time.time()
            print("Generating private key using GGH's matrix transformations tecnique...")
            l = n
            k = fmpz(l * math.ceil(math.sqrt(self.dimension) + 1))
            
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
        else:
            time_start = time.time()
            print("Generating private key using Micciancio's random matrices tecnique...")
            time_start = time.time()
            while True:
                try:
                    R = self.numpy_to_fmpz_mat(np.random.randint(-n, n, size=(n,n)))
                    R_inv = R.inv()
                except Exception as e:
                    print(e)
                    continue
                else:
                    break
            R = R.lll()
            R_inv = R.inv()
        
        if self.debug:
            print(f"..time taken: {time.time() - time_start}")
        self.good_basis = R
        print("Generating public key...")

        self.bad_basis = H = R.hnf()

        self.public_key = H

        self.private_key = (R_inv, R)

    def encrypt(self):

        if self.x is None:
            x = self.reduce_mod_B()
        
        H = self.bad_basis
        r = self.error

        self.ciphertext = r - x * H

    def decrypt(self):
        R_inv, R = self.private_key
        
        c = self.ciphertext
        x = c * R_inv

        rounded_x = fmpz_mat(x.nrows(), x.ncols())

        for i in range(x.nrows()):
            for j in range(x.ncols()):
                rounded_x[i, j] = round(x[i,j])
                                
        result = rounded_x * R
        
        return c - result
    
    def get_hadamard_ratio(self, basis = None):
        matrix = basis
        norms = []
        
        for i in range(matrix.nrows()):
            row = [matrix[i, j] for j in range(matrix.ncols())]
            norm = self.vector_norm(row)
            norms.append(norm)
           
        denominator = math.prod(norms)
        numerator = abs(Decimal(matrix.det().str()))

        result = (numerator / denominator) ** Decimal(1 / self.dimension)
        return result
    
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