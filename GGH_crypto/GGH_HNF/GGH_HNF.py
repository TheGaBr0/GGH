import sympy as sp
import random
import math
from flint import fmpz_mat, fmpq_mat, fmpq, fmpz
from decimal import Decimal, getcontext
import time

from ..Utils.Utils import Utils



class GGHHNFCryptosystem:
    def __init__(self, dimension, R=None, B=None, x=None, e=None, rho_check=True, random_private=True, debug=True):
        self.dimension = dimension
        
        self.x = x
        self.ciphertext = None
        self.error = e

        self.private_basis = R
        self.public_basis = B
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
                raise ValueError(f"[GGH-HNF] Private basis must be a {dimension}x{dimension} matrix, but got {R.nrows()}")
            
            self.generate_keys_from_R_or_B()
        else:
            self.generate_keys()
            
        if e is None:
            self.generate_random_error()
        else:
            if self.debug:
                print(f"[GGH-HNF] Length of error vector is: {Utils.vector_l2_norm(e)}")

    def generate_keys_from_R_or_B(self):

        if self.debug:
            print("[GGH-HNF] Private basis given as input, inverting it..")
    
        R = self.private_basis
        R_inv = R.inv()

        if self.rho_check:
            if self.debug:
                print("[GGH-HNF] Calculating rho...")
            self.R_rho = self.calculate_rho(self.private_basis)

        if self.public_basis is None:
            if self.debug:
                print("[GGH-HNF] Generating public basis using HNF of the good basis...")
            H = R.hnf()
            self.public_basis = H
        else:
            if self.debug:
                print("[GGH-HNF] Using the provided public basis as the public key")
            H = self.public_basis
        
        self.public_key = (H, self.R_rho)
        self.private_key = (R_inv, R)

    def min_norm_row(self, matrix):
        norms = []
        for j in range(matrix.nrows()):
            row = [matrix[j, i] for i in range(matrix.ncols())]
            norms.append(Utils.vector_l2_norm(row))
        min_norm = min(norms)
        return norms[norms.index(min_norm)]

    def calculate_rho(self, basis):
        if self.debug:
            print("[GGH-HNF] Gram Schmidt orthogonalization started...")
        time_start = time.time()
        basis_orthogonalized = Utils.gram_schmidt(basis)
        if self.debug:
            print(f"[GGH-HNF] Time taken: {time.time() - time_start}")
        min_norm = self.min_norm_row(basis_orthogonalized)
        rho = Decimal(0.5) * min_norm

        if self.debug:
            print(f'[GGH-HNF] Rho is: {rho:f}')
        return rho
    
    def generate_random_error(self):
        n = self.dimension
        if self.rho_check:
            rand_val = n
            random_elements = [[random.randint(-rand_val, rand_val) for _ in range(n)]]
            error = fmpz_mat(random_elements)
            error_norm = Utils.vector_l2_norm(error)
            while(not error_norm < Decimal(0.9) * self.R_rho):
                rand_val -= 1
                random_elements = [[random.randint(-rand_val, rand_val) for _ in range(n)]]
                error = fmpz_mat(random_elements)
                error_norm = Utils.vector_l2_norm(error)
            self.error = error
        else:
            random_elements = [[random.randint(-28, 28) for _ in range(n)]]
            self.error = fmpz_mat(random_elements)
            error_norm = Utils.vector_l2_norm(self.error)
        
        if self.debug:
            print(f"[GGH-HNF] Length of error vector is: {error_norm}")
    
    def reduce_mod_B(self):
        r = fmpq_mat(self.error)
        H = self.public_basis
        x = fmpz_mat(r.nrows(), r.ncols())
        
        # Iterate over each coordinate starting from the last one
        for i in reversed(range(r.ncols())):
            sum_j = sum(H[j, i] * x[0, j] for j in range(i + 1, self.dimension))
            # Compute x[i] using the formula
            x[0, i] = math.floor((r[0, i] - sum_j) / H[i, i])

        return x

    def generate_keys(self):
        n = self.dimension
        tries = 0
        if not self.random_private:
            if self.debug:
                print("[GGH-HNF] Generating private basis using GGH's matrix transformations technique...")
            time_start = time.time()
            l = n
            k = fmpz(l * math.ceil(math.sqrt(self.dimension) + 1))
            
            while True:
                try:
                    R = fmpz_mat([[random.randint(-l, l-1) for _ in range(self.dimension)] for _ in range(self.dimension)])
                    I = Utils.sympy_to_fmpz_mat(sp.eye(self.dimension))
                    KI = k * I
                    R += KI
                    R_inv = R.inv()
                    tries += 1
                except Exception as e:
                    print(e)
                    continue
                else:
                    break
        else:
            if self.debug:
                print("[GGH-HNF] Generating private basis using Micciancio's random matrices tecnique...")
            time_start = time.time()
            while True:
                try:
                    R = fmpz_mat([[random.randint(-n, n - 1) for _ in range(n)] for _ in range(n)])
                    R_inv = R.inv()
                    tries += 1
                except Exception as e:
                    print(e)
                    continue
                else:
                    break
            R = R.lll()
            R_inv = R.inv()
        
        if self.debug:
            print(f"[GGH-HNF] Time taken: {time.time() - time_start} with {tries} tries")
        self.private_basis = R

        if self.rho_check:
            if self.debug:
                print("[GGH-HNF] Calculating rho...")
            self.R_rho = self.calculate_rho(self.private_basis)

        if self.debug:
            print("[GGH-HNF] Generating public basis...")
        time_start = time.time()
        self.public_basis = H = R.hnf()
        if self.debug:
            print(f"[GGH-HNF] Time taken: {time.time() - time_start}")

        self.public_key = (H, self.R_rho)

        self.private_key = (R_inv, R)

    def encrypt(self):
        if self.debug:
            print(f"[GGH-HNF] Encrypting...")
        time_start = time.time()
        if self.x is None:
            x = self.reduce_mod_B()
        
        H = self.public_basis
        r = self.error

        self.ciphertext = r - x * H
        if self.debug:
            print(f"[GGH] Time taken: {time.time() - time_start}")

    def decrypt(self):
        if self.debug:
            print(f"[GGH] Decrypting...")
        time_start = time.time()
        R_inv, R = self.private_key
        
        c = self.ciphertext
        x = c * R_inv

        rounded_x = fmpz_mat(x.nrows(), x.ncols())

        for i in range(x.nrows()):
            for j in range(x.ncols()):
                rounded_x[i, j] = round(x[i,j])
                                
        result = rounded_x * R
        
        if self.debug:
            print(f"[GGH] Time taken: {time.time() - time_start}")


        return c - result
    
    
