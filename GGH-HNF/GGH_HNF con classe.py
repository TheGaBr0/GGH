from GGH_HNF import GGHHNFCryptosystem
import sympy as sp
import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpq

dimension = 100


R = fmpz_mat([[12, -4, -1], [1,  8, -1], [-4,  1, 14]])
e = fmpz_mat([[1,1,2]])

GGH_object = GGHHNFCryptosystem(dimension = dimension, rho_check=True, random_private=False)
GGH_object.encrypt()

B = GGH_object.public_key

R = GGH_object.private_key[1]
    
error = GGH_object.error

ciphertext = GGH_object.ciphertext

decrypted_message = GGH_object.decrypt()

print(f"error: {error}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"error: {decrypted_message}")




