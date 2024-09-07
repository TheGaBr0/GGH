import sympy as sp
import random
import math
from flint import fmpz_mat, fmpq_mat, fmpq, fmpz
from decimal import Decimal, getcontext
import time
import numpy as np

from ..Utils.Utils import Utils

class GGHHNFCryptosystem:
    def __init__(self, dimension, private_basis=None, public_basis=None, lattice_point=None, error=None, GGH_private=False, debug=False):
        self.dimension = dimension
        
        self.lattice_point = lattice_point
        self.ciphertext = None
        self.error = error

        self.private_basis = private_basis
        self.public_basis = public_basis
        self.unimodular = None
        
        self.R_rho = None
        self.GGH_private = GGH_private
        self.debug = debug

        self.private_key = None
        self.public_key = None

        if private_basis is not None:
            if private_basis.nrows() != dimension or private_basis.ncols() != dimension:
                raise ValueError(f"[GGH-HNF] Private basis must be a {dimension}x{dimension} matrix, but got {private_basis.nrows()}x{private_basis.ncols()}")
        
        if public_basis is not None:
            if public_basis.nrows() != dimension or public_basis.ncols() != dimension:
                raise ValueError(f"[GGH-HNF] Public basis must be a {dimension}x{dimension} matrix, but got {public_basis.nrows()}x{public_basis.ncols()}")

        if lattice_point is not None:
            if lattice_point.ncols() != dimension:
                raise ValueError(f"[GGH-HNF] Lattice point vector must have length {dimension}, but got length {lattice_point.ncols()}")
        
        if error is not None:
            if error.ncols() != dimension:
                raise ValueError(f"[GGH-HNF] Error vector must have length {dimension}, but got length {error.ncols()}")

        if self.private_basis is not None:
            self.generate_keys_from_R_or_B()
        else:
            self.generate_keys()
            
        if self.error is None:
            self.generate_random_error()
        else:
            if self.debug:
                print(f"[GGH-HNF] Length of error vector is: {Utils.vector_l2_norm(self.error)}")

    def generate_keys_from_R_or_B(self):

        if self.debug:
            print("[GGH-HNF] Private basis given as input, inverting it..")
    
        R = self.private_basis
        R_inv = R.inv()

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
        basis_orthogonalized = Utils.gram_schmidt(np.array(basis.tolist()).astype(int))
        basis_orthogonalized = Utils.npsp_to_fmpq_mat(basis_orthogonalized)
        min_norm = self.min_norm_row(basis_orthogonalized)
        rho = Decimal(0.5) * min_norm

        if self.debug:
            print(f'[GGH-HNF] Rho is: {rho:f}')
        return rho
    
    def generate_random_error(self):
        n = self.dimension
        max_norm = Decimal(0.75) * self.R_rho

        while True:
            error = fmpz_mat([[random.randint(-n, n) for _ in range(self.dimension)]])
      
            error_norm = Utils.vector_l2_norm(error)
            
            if error_norm < max_norm:
                break
            n -= 1
        
        self.error = error

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
        if self.GGH_private:
            if self.debug:
                print("[GGH-HNF] Generating private basis using GGH's matrix transformations technique...")
                time_start = time.time()
            l = n
            k = fmpz(l * math.ceil(math.sqrt(self.dimension) + 1))
            
            while True:
                try:
                    R = fmpz_mat([[random.randint(-l, l-1) for _ in range(self.dimension)] for _ in range(self.dimension)])
                    I = Utils.npsp_to_fmpz_mat(sp.eye(self.dimension))
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
                    R = R.lll()
                    R_inv = R.inv()
                    tries += 1
                except Exception as e:
                    print(e)
                    continue
                else:
                    break
        
        if self.debug:
            priv_time = time.time() - time_start
            print(f"[GGH-HNF] Time taken: {priv_time} with {tries} tries")
        self.private_basis = R

        if self.debug:
            time_start = time.time()
            print("[GGH-HNF] Calculating rho...")
        self.R_rho = self.calculate_rho(self.private_basis)
        if self.debug:
            rho_time = time.time() - time_start
            print(f"[GGH-HNF] Time taken: {rho_time}")

        if self.debug:
            print("[GGH-HNF] Generating public basis...")
            time_start = time.time()
        self.public_basis = H = R.hnf()
        if self.debug:
            pub_time = time.time() - time_start
            print(f"[GGH-HNF] Time taken: {pub_time}")

        self.public_key = (H, self.R_rho)
        self.private_key = R

    def encrypt(self):
        if self.debug:
            print(f"[GGH-HNF] Encrypting...")
            time_start = time.time()
        if self.lattice_point is None:
            x = self.reduce_mod_B()
        
        H = self.public_basis
        r = self.error

        self.ciphertext = r - x * H
        if self.debug:
            enc_time = time.time() - time_start
            print(f"[GGH-HNF] Time taken: {enc_time}")

    def decrypt(self):
        if self.debug:
            print(f"[GGH-HNF] Decrypting...")
            time_start = time.time()
        
        CVP = Utils.babai_rounding(self.private_basis, self.ciphertext)

        if self.debug:
            dec_time = time.time() - time_start
            print(f"[GGH-HNF] Time taken: {dec_time}")


        return fmpq_mat(self.ciphertext) - CVP
    
    
