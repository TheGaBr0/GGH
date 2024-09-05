import random
import sympy as sp
import math
from flint import fmpz_mat, fmpz, fmpq, fmpq_mat
from decimal import Decimal, getcontext
from fractions import Fraction
import time

from ..Utils.Utils import Utils

class GGHCryptosystem:
    def __init__(self, dimension, private_basis=None, public_basis=None, unimodular=None, message=None, 
                 error=None, sigma = None, integer_sigma=True, debug=False):
        self.dimension = dimension
        self.integer_sigma = integer_sigma
        self.sigma = sigma
        
        self.message = message
        self.ciphertext = None
        self.error = error

        self.private_basis = private_basis
        self.public_basis = public_basis
        self.unimodular = unimodular
        self.debug = debug

        self.private_key = None
        self.public_key = None

        # Check dimensions of bases if provided
        if private_basis is not None:
            if private_basis.nrows() != dimension or private_basis.ncols() != dimension:
                raise ValueError(f"[GGH] Private basis must be a {dimension}x{dimension} matrix, but got {private_basis.nrows()}x{private_basis.ncols()}")
        
        if public_basis is not None:
            if public_basis.nrows() != dimension or public_basis.ncols() != dimension:
                raise ValueError(f"[GGH] Public basis must be a {dimension}x{dimension} matrix, but got {public_basis.nrows()}x{public_basis.ncols()}")
        
        if unimodular is not None:
            if unimodular.nrows() != dimension or unimodular.ncols() != dimension:
                raise ValueError(f"[GGH] Unimodular matrix must be a {dimension}x{dimension} matrix, but got {unimodular.nrows()}x{unimodular.ncols()}")

        # Check dimensions of vectors if provided
        if message is not None:
            if message.ncols() != dimension:
                raise ValueError(f"[GGH] Message vector must have length {dimension}, but got length {message.ncols()}")
        
        if error is not None:
            if error.ncols() != dimension:
                raise ValueError(f"[GGH] Error vector must have length {dimension}, but got length {error.ncols()}")


        if private_basis is not None:
            self.generate_keys_from_R_or_B()
        else:
            self.generate_keys()

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

        return Utils.sympy_to_fmpz_mat(T)

    def generate_sigma(self, R_inv):
        getcontext().prec = 50
        rho = Utils.vector_l1_norm(R_inv)

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
            random_elements = [fmpq(Fraction(item).numerator, Fraction(item).denominator) for item in random_elements]
            self.error = fmpq_mat([random_elements])

    def generate_random_message(self):
        random_elements = [random.randint(-128, 127) for _ in range(self.dimension)]
        self.message = fmpq_mat([random_elements])

    def generate_keys_from_R_or_B(self):

        if self.debug:
            print("[GGH] Private basis given as input, inverting it..")

        R_inv = self.private_basis.inv()

        if not self.sigma:
            self.sigma = self.generate_sigma(R_inv)

        if self.public_basis is None:
            if self.debug:
                print("[GGH] Generating public key...")

            if not self.unimodular:
                self.unimodular = self.random_unimodular(self.dimension, 2)
            self.public_basis = self.unimodular * self.private_basis
        else:
            if self.debug:
                print("[GGH] Using the provided public basis as the public key")
        
        self.public_key = (self.public_basis, self.sigma)
        self.private_key = (R_inv, self.private_basis)

    def generate_keys(self):
        tries = 0
        l = 4
        k = fmpz(l * math.ceil(math.sqrt(self.dimension) + 1))
        if self.debug:
            print("[GGH] Generating private basis...")
            time_start = time.time()
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

        if self.debug:
            priv_time = time.time() - time_start
            print(f"[GGH] Time taken: {priv_time} with {tries} tries")

        self.private_basis = R

        if not self.sigma:
            if self.debug:
                print("[GGH] Generating sigma...")
                time_start = time.time()
            self.sigma = self.generate_sigma(R_inv)
            if self.debug:
                sigma_time = time.time() - time_start
                print(f"Generated sigma is {self.sigma}, time taken: {sigma_time}")

        if not self.unimodular:
            if self.debug:
                print(f"Generating unimodular matrix...")
                time_start = time.time()
            self.unimodular = self.random_unimodular(self.dimension, 2)
            if self.debug:
                unim_time = time.time() - time_start
                print(f"[GGH] Time taken: {unim_time}")


        if self.debug:
            print(f"[GGH] Generating public basis...")
            time_start = time.time()
        self.public_basis = self.unimodular * self.private_basis
        if self.debug:
            pub_time = time.time() - time_start
            print(f"[GGH] Time taken: {pub_time}")

        self.public_key = (self.public_basis, self.sigma)
        self.private_key = (R_inv, R)

    def encrypt(self):
        if self.debug:
            print(f"[GGH] Encrypting...")
            time_start = time.time()
        if self.error is None:
            self.generate_error()

        if self.message is None:
            self.generate_random_message()

        B = self.public_key[0]

        self.ciphertext = self.message * B + self.error
        if self.debug:
            enc_time = time.time() - time_start
            print(f"[GGH] Time taken: {enc_time}")

    def decrypt(self):
        if self.debug:
            print(f"[GGH] Decrypting...")
            time_start = time.time()
        
        CVP = Utils.babai_rounding(self.private_basis, self.ciphertext)
        
        result = CVP * self.public_basis.inv()

        if self.debug:
            dec_time = time.time() - time_start
            print(f"[GGH] Time taken: {dec_time}")

        return result
