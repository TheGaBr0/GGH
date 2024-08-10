import numpy as np
import random
import sympy as sp
import math
from flint import fmpz_mat, fmpz, fmpq, fmpq_mat
from decimal import Decimal, getcontext

class GGHCryptosystem:
    def __init__(self, dimension, R=None, B=None, m=None, e=None, integer_sigma=True):
        self.dimension = dimension
        self.integer_sigma = integer_sigma
        
        self.message = m
        self.ciphertext = None
        self.error = e

        self.private_basis = R
        self.public_basis = B
        self.unimodular = None

        self.private_key = None
        self.public_key = None

        if R is not None:
            # Check if R is a 2D array or matrix
            if R.nrows() != dimension:
                raise ValueError(f"R must be a {dimension}x{dimension} matrix, but got {R.nrows()}")
            
            self.private_basis = R
            if B is None:
                self.generate_keys_from_R()
            else:
                self.public_basis = B
                self.generate_keys_with_B()
        else:
            self.generate_keys()

        if m is None:
            self.generate_random_message()

        sp.init_printing(use_unicode=True)
        

    def sympy_to_fmpz_mat(self, basis_sympy):
        return fmpz_mat([[int(item) for item in sublist] for sublist in basis_sympy.tolist()])
    
    def numpy_to_fmpz_mat(self, numpy_matrix):
        return fmpz_mat([[int(item) for item in sublist] for sublist in numpy_matrix.tolist()])

    def random_unimodular(self, dim, mix):
        T = sp.eye(dim)
        while mix > 0:
            rows = list(range(dim))
            random.shuffle(rows)
            for i in rows:
                x = sp.zeros(1, dim)
                x[0, i] = 1
                for k in range(dim):
                    if k != i:
                        x[0, k] = np.random.choice([-1, 0, 1], p=[1/7, 5/7, 1/7])
                T[i, :] = x * T
            mix -= 1
        return self.sympy_to_fmpz_mat(T)

    def generate_sigma(self, R_inv):
        rho = fmpq(max(sum(abs(x) for x in row) for row in R_inv.tolist())) #L1 norm

        sigma_max = 1 / (2 * rho)

        if self.integer_sigma:
            return int(math.floor(sigma_max)) if math.floor(sigma_max) < 3 else 3 #standard sigma is 3
        else:
            return sigma_max


    def generate_error(self):
        sigma = self.public_key[1]

        random_elements = [random.choice([-sigma, sigma]) for _ in range(self.dimension)]
        
        if isinstance(sigma, int):
            self.error = fmpz_mat([random_elements])
        else:
            random_elements = [fmpq(sp.Rational(item).numerator, sp.Rational(item).denominator) for item in random_elements]
            self.error = fmpq_mat([random_elements])

    def generate_random_message(self):
        random_elements = [random.randint(-128, 127) for _ in range(self.dimension)]
        self.message = fmpq_mat([random_elements])

    def generate_keys_from_R(self):
        R = self.private_basis
        R_inv = R.inv()
        sigma = self.generate_sigma(R_inv)
        print("Generating public key...")
        U = self.random_unimodular(self.dimension, 12)
        self.unimodular = U
        self.public_basis = B = U * R
        print("Computing results...")
        self.public_key = (B, sigma)
        self.private_key = (R_inv, R)

    def generate_keys_with_B(self):
        R = self.private_basis
        R_inv = R.inv()
        sigma = self.generate_sigma(R_inv)
        print("Computing results...")
        self.public_key = (self.public_basis, sigma)
        self.private_key = (R_inv, R)

    def generate_keys(self):
        l = 4
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
        self.private_basis = R
        sigma = self.generate_sigma(R_inv)
        print("Generating public key...")
        U = self.random_unimodular(self.dimension, 12)
        self.unimodular = U
        self.public_basis = B = U * R
        print("Computing results...")
        self.public_key = (B, sigma)
        self.private_key = (R_inv, R)

    def encrypt(self):
        if self.error is None:
            self.generate_error()
        B = self.public_key[0]
        self.ciphertext = self.message * B + self.error

    def decrypt(self):
        R_inv, R = self.private_key
        B_inv = self.public_basis.inv()
        c = self.ciphertext
        a = c * R_inv
        
        cols = a.ncols()
        rows = a.nrows()

        for i in range(rows):
            for j in range(cols):
                a[i,j] = round(a[i,j])
        
        print(a)

        result = a * R * B_inv
        return result
    
    def row_norm(self, row):
        return Decimal(sum(Decimal(int(x))**2 for x in row)).sqrt()

    def get_hadamard_ratio(self, basis = None):
        
        matrix = basis
        norms = []
        
        for i in range(matrix.nrows()):
            row = [matrix[i, j] for j in range(matrix.ncols())]
            norm = self.row_norm(row)
            norms.append(norm)
           
        
        denominator = math.prod(norms)
        numerator = abs(Decimal(matrix.det().str()))


        result = (numerator / denominator) ** Decimal(1 / self.dimension)
        return result