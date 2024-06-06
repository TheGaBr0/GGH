import numpy as np
import random
import sympy as sp
import math
from flint import fmpz_mat, fmpz, fmpq, fmpq_mat
from decimal import Decimal, getcontext

class GGHCryptosystem:
    def __init__(self, dimension, use_random=True, integer_sigma=True):
        self.dimension = dimension
        self.integer_sigma = integer_sigma
        
        self.message = None
        self.ciphertext = None
        self.error = None

        self.private_basis = None
        self.public_basis = None
        self.unimodular = None

        self.private_key = None
        self.public_key = None

        if use_random:
            self.generate_random_message()
            self.generate_keys()

        sp.init_printing(use_unicode=True)
        

    def sympy_to_fmpz_mat(self, basis_sympy):
        return fmpz_mat([[int(item) for item in sublist] for sublist in basis_sympy.tolist()])
    
    def numpy_to_fmpz_mat(self, numpy_matrix):
        return fmpz_mat([[int(item) for item in sublist] for sublist in numpy_matrix.tolist()])

    def random_unimodular(self, dim, mix):
        T = sp.eye(dim)
        while mix > 0:
            cols = list(range(dim))
            random.shuffle(cols)
            for j in cols:
                x = sp.zeros(dim, 1)
                x[j, 0] = 1
                for k in range(dim):
                    if k != j:
                        x[k, 0] = np.random.choice([-1, 0, 1], p=[1/7, 5/7, 1/7])
                T[:, j] = T * x
            mix -= 1
        return self.sympy_to_fmpz_mat(T)

    def generate_sigma(self, R_inv):

        #rho = fmpq(max(sum(abs(x) for x in row) for row in R_inv.tolist())) #infinite norm
                
        rho = fmpq(max(sum(abs(x) for x in row) for row in R_inv.transpose().tolist())) #L1 norm, i need to transpose the matrix cause i want the absolute column sums

        sigma_max = 1 / (2 * rho)

        if self.integer_sigma:
            return int(math.floor(sigma_max)) if math.floor(sigma_max) < 3 else 3 #standard sigma is 3
        else:
            return sigma_max


    def generate_error(self):
        sigma = self.public_key[1]

        random_elements = [random.choice([-sigma, sigma]) for _ in range(self.dimension)]
        
        if isinstance(sigma, int):
            self.error = fmpz_mat([random_elements]).transpose()
        else:
            random_elements = [fmpq(sp.Rational(item).numerator, sp.Rational(item).denominator) for item in random_elements]
            self.error = fmpq_mat(random_elements).transpose()

    def generate_random_message(self):
        random_elements = [random.randint(-128, 127) for _ in range(self.dimension)]
        self.message = fmpq_mat([random_elements]).transpose()

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
        U = self.random_unimodular(self.dimension, 2)
        self.unimodular = U
        self.public_basis = B = R * U
        print("Computing results...")
        self.public_key = (B, sigma)
        B_inv = B.inv()
        self.private_key = (R_inv, B_inv * R)

    def encrypt(self):
        self.generate_error()
        B = self.public_key[0]
        self.ciphertext = B * self.message + self.error

    def decrypt(self):
        R_inv, U_inv = self.private_key
        c = self.ciphertext
        a = R_inv * c
        
        cols = a.ncols()
        rows = a.nrows()

        for i in range(rows):
            for j in range(cols):
                a[i,j] = round(a[i,j])
                
        result = U_inv * a
        return result
    
    def column_norm(self, col):
        return Decimal(sum(Decimal(int(x))**2 for x in col)).sqrt()

    def get_hadamard_ratio(self, basis = None):
        
        matrix = basis
        norms = []
        
        for j in range(matrix.ncols()):
            column = [matrix[i, j] for i in range(matrix.nrows())]
            norm = self.column_norm(column)
            norms.append(norm)
           
        
        denominator = math.prod(norms)
        numerator = abs(Decimal(matrix.det().str()))


        result = (numerator / denominator) ** Decimal(1 / self.dimension)
        return result
