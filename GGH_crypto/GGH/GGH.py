import random
import sympy as sp
import math
from flint import fmpz_mat, fmpz, fmpq, fmpq_mat
from decimal import Decimal, getcontext
import time

class GGHCryptosystem:
    def __init__(self, dimension, R=None, B=None, m=None, e=None, integer_sigma=True, debug=True):
        self.dimension = dimension
        self.integer_sigma = integer_sigma
        
        self.message = m
        self.ciphertext = None
        self.error = e

        self.private_basis = R
        self.public_basis = B
        self.unimodular = None
        self.debug = debug

        self.private_key = None
        self.public_key = None

        if R is not None:
            # Check if R is a 2D array or matrix
            if R.nrows() != dimension:
                raise ValueError(f"[GGH] Private basis must be a {dimension}x{dimension} matrix, but got {R.nrows()}")
            
            self.generate_keys_from_R_or_B()
        else:
            self.generate_keys()

        if m is None:
            self.generate_random_message()

        sp.init_printing(use_unicode=True)
        

    def sympy_to_fmpz_mat(self, basis_sympy):
        return fmpz_mat([[int(item) for item in sublist] for sublist in basis_sympy.tolist()])

    def random_unimodular(self, dim, mix):
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

        return self.sympy_to_fmpz_mat(T)

    def generate_sigma(self, R_inv):
        getcontext().prec = 50
        rho = max(sum(abs(Decimal(int(x.numer())) / Decimal(int(x.denom()))) for x in row) for row in R_inv.tolist()) #L1 norm

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

    def generate_keys_from_R_or_B(self):

        if self.debug:
            print("[GGH] Private basis given as input, inverting it..")
    
        R = self.private_basis
        R_inv = R.inv()
        sigma = self.generate_sigma(R_inv)

        if self.public_basis is None:
            if self.debug:
                print("[GGH] Generating public key...")
            U = self.random_unimodular(self.dimension, 2)
            self.unimodular = U
            self.public_basis = B = U * R
        else:
            if self.debug:
                print("[GGH] Using the provided public basis as the public key")
            B = self.public_basis
        
        self.public_key = (B, sigma)
        self.private_key = (R_inv, R)

    def generate_keys(self):
        if self.debug:
            print("[GGH] Generating private basis...")
        tries = 0
        time_start = time.time()
        l = 4
        k = fmpz(l * math.ceil(math.sqrt(self.dimension) + 1))
        
        while True:
            try:
                R = fmpz_mat([[random.randint(-l, l-1) for _ in range(self.dimension)] for _ in range(self.dimension)])
                I = self.sympy_to_fmpz_mat(sp.eye(self.dimension))
                KI = k * I
                R += KI
                R_inv = R.inv()
                tries += 1
            except Exception as e:
                print(e)
                continue
            else:
                break

        if self.debug:
            print(f"[GGH] Time taken: {time.time() - time_start} with {tries} tries")

        self.private_basis = R
        sigma = self.generate_sigma(R_inv)

        if self.debug:
            print(f"Generated sigma is {sigma}")

        if self.debug:
            print(f"Generating unimodular matrix")
        time_start = time.time()
        U = self.random_unimodular(self.dimension, 2)

        if self.debug:
            print(f"[GGH] Time taken: {time.time() - time_start}")

        self.unimodular = U
        if self.debug:
            print(f"[GGH] Generating public basis")
        time_start = time.time()
        self.public_basis = B = U * R
        if self.debug:
            print(f"[GGH] Time taken: {time.time() - time_start}")

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
        
        result = a * R * B_inv
        return result
    
    def vector_norm(self, row):
        if isinstance(row, fmpz_mat):
            row = fmpq_mat(row)
        getcontext().prec = 50
        return Decimal(sum((Decimal(int(x.numer())) / Decimal(int(x.denom()))) ** 2 for x in row)).sqrt()

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