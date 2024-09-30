import sympy as sp
import random
import math
from flint import fmpz_mat, fmpq_mat, fmpq, fmpz
from decimal import Decimal
import time
import numpy as np
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

from ..Utils.Utils import Utils

class GGHHNFCryptosystem:
    """
    Implementation of the GGH (Goldreich-Goldwasser-Halevi) cryptographic system with Hermite Normal Form (HNF) modifications.

    This class provides methods for generating keys, encrypting, and decrypting messages
    using the GGH lattice-based cryptographic system with HNF optimizations.

    Attributes:
        dimension (int): The dimension of the lattice.
        lattice_point (fmpz_mat): The lattice point to be encrypted.
        ciphertext (fmpq_mat): The encrypted text.
        error (fmpz_mat): The error vector.
        alpha (float): The alpha parameter for error generation.
        private_basis (fmpz_mat): The private basis of the lattice.
        public_basis (fmpz_mat): The public basis of the lattice.
        R_rho (Decimal): The rho parameter of the private basis.
        GGH_private (bool): If True, uses GGH's matrix transformation technique for private basis generation.
        debug (bool): If True, prints debug information.
        private_key (tuple): The private key (R_inv, R).
        public_key (tuple): The public key (H, R_rho).

    Args:
        dimension (int): The dimension of the lattice.
        private_basis (fmpz_mat, optional): A predefined private basis.
        public_basis (fmpz_mat, optional): A predefined public basis.
        lattice_point (fmpz_mat, optional): A predefined lattice point.
        error (fmpz_mat, optional): A predefined error vector.
        alpha (float, optional): The alpha parameter for error generation. Default is 0.75.
        GGH_private (bool, optional): Whether to use GGH's matrix transformation technique. Default is False.
        debug (bool, optional): Whether to enable debug prints. Default is False.

    Raises:
        ValueError: If the dimensions of the provided bases or vectors do not match.
    """
    def __init__(self, dimension, private_basis=None, public_basis=None, lattice_point=None, error=None, alpha=0.75, GGH_private=False, debug=False):
        self.dimension = dimension
        
        self.lattice_point = lattice_point
        self.ciphertext = None
        self.error = error
        self.alpha = alpha

        self.private_basis = private_basis
        self.public_basis = public_basis
        
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
                logger.info(f"[GGH-HNF] Length of error vector is: {Utils.vector_l2_norm(self.error)}")

    def generate_keys_from_R_or_B(self):
        """
        Generates keys from a provided private or public basis.
        """
        if self.debug:
            logger.info("[GGH-HNF] Private basis given as input, inverting it..")
    
        R = self.private_basis
        R_inv = R.inv()

        if self.debug:
            logger.info("[GGH-HNF] Calculating rho...")
            self.R_rho = self.calculate_rho(self.private_basis)
            

        if self.public_basis is None:
            if self.debug:
                logger.info("[GGH-HNF] Generating public basis using HNF of the good basis...")
            H = R.hnf()
            self.public_basis = H
        else:
            if self.debug:
                logger.info("[GGH-HNF] Using the provided public basis as the public key")
            H = self.public_basis
        
        self.public_key = (H, self.R_rho)
        self.private_key = (R_inv, R)

    def min_norm_row(self, matrix):
        """
        Finds the minimum norm row in the given matrix.

        Args:
            matrix (fmpz_mat or fmpq_mat): The input matrix.

        Returns:
            float: The minimum norm among all rows.
        """
        norms = []
        for j in range(matrix.nrows()):
            row = [matrix[j, i] for i in range(matrix.ncols())]
            norms.append(Utils.vector_l2_norm(row))
        min_norm = min(norms)
        return norms[norms.index(min_norm)]

    def calculate_rho(self, basis):
        """
        Calculates the rho parameter for the given basis.

        Args:
            basis (fmpz_mat): The input basis.

        Returns:
            Decimal: The calculated rho value.
        """
        basis_orthogonalized = Utils.gram_schmidt(np.array(basis.tolist()).astype(int))
        basis_orthogonalized = Utils.npsp_to_fmpq_mat(basis_orthogonalized)
        min_norm = self.min_norm_row(basis_orthogonalized)
        rho = Decimal(0.5) * min_norm

        if self.debug:
            logger.info(f'[GGH-HNF] Rho is: {rho:f}')
        return rho
    
    def generate_random_error(self):
        """
        Generates a random error vector based on the alpha and R_rho parameters.
        """
        n = self.dimension
        max_norm = Decimal(self.alpha) * self.R_rho

        while True:
            error = fmpz_mat([[random.randint(-n, n) for _ in range(self.dimension)]])
      
            error_norm = Utils.vector_l2_norm(error)
            
            if error_norm < max_norm:
                break
            n -= 1
        
        self.error = error

        if self.debug:
            logger.info(f"[GGH-HNF] Length of error vector is: {error_norm}")
    
    def reduce_mod_B(self):
        """
        Generates the private and public keys of the cryptographic system.
        """
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
        """
        Encrypts the lattice point using the public key and an error vector.
        """
        n = self.dimension
        tries = 0
        if self.GGH_private:
            if self.debug:
                logger.info("[GGH-HNF] Generating private basis using GGH's matrix transformations technique...")
                time_start = time.time()
            l = n
            k = fmpz(l * math.ceil(math.sqrt(self.dimension) + 1))
            
            while True:
                R = fmpz_mat([[random.randint(-l, l-1) for _ in range(self.dimension)] for _ in range(self.dimension)])
                I = Utils.npsp_to_fmpz_mat(sp.eye(self.dimension))
                KI = k * I
                R += KI
                
                tries += 1

                if R.det() != 0:
                    break
        else:
            if self.debug:
                logger.info("[GGH-HNF] Generating private basis using Micciancio's random matrices tecnique...")
                time_start = time.time()
            while True:
                R = fmpz_mat([[random.randint(-n, n - 1) for _ in range(n)] for _ in range(n)])
                R = R.lll()
                
                tries += 1

                if R.det() != 0:
                    break
               
        
        if self.debug:
            priv_time = time.time() - time_start
            logger.info(f"[GGH-HNF] Time taken: {priv_time} with {tries} tries")
        self.private_basis = R

        if self.debug:
            time_start = time.time()
            logger.info("[GGH-HNF] Calculating rho...")
        self.R_rho = self.calculate_rho(self.private_basis)
        if self.debug:
            rho_time = time.time() - time_start
            logger.info(f"[GGH-HNF] Time taken: {rho_time}")

        if self.debug:
            logger.info("[GGH-HNF] Generating public basis...")
            time_start = time.time()
        self.public_basis = H = R.hnf()
        if self.debug:
            pub_time = time.time() - time_start
            logger.info(f"[GGH-HNF] Time taken: {pub_time}")

        self.public_key = (H, self.R_rho)
        self.private_key = R

    def encrypt(self):
        if self.debug:
            logger.info(f"[GGH-HNF] Encrypting...")
            time_start = time.time()
        if self.lattice_point is None:
            x = self.reduce_mod_B()
        
        H = self.public_basis
        r = self.error

        self.ciphertext = r - x * H
        if self.debug:
            enc_time = time.time() - time_start
            logger.info(f"[GGH-HNF] Time taken: {enc_time}")

    def decrypt(self):
        """
        Decrypts the ciphertext using the private key.

        Returns:
            fmpq_mat: The decrypted lattice point.
        """
        if self.debug:
            logger.info(f"[GGH-HNF] Decrypting...")
            time_start = time.time()
        
        CVP = Utils.babai_rounding(self.private_basis, self.ciphertext)

        if self.debug:
            dec_time = time.time() - time_start
            logger.info(f"[GGH-HNF] Time taken: {dec_time}")


        return fmpq_mat(self.ciphertext) - CVP
    
    