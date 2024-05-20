from GGH_HNF import GGHCryptosystem
import sympy as sp
import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpz
import os
import subprocess
import symengine as sym
import time
dimension = 400
GGH_object = GGHCryptosystem(dimension = dimension)
GGH_object.encrypt()

R = GGH_object.private_key[1]
H = GGH_object.public_key
    
print(GGH_object.generate_rho(R))

print(GGH_object.get_hadamard_ratio())

print(GGH_object.generate_rho(H))

print(GGH_object.generate_rho(GGH_object.error))

H = GGH_object.public_key
    
error = GGH_object.error

ciphertext = GGH_object.ciphertext

decrypted_message = GGH_object.decrypt()

print(f"error: {error.transpose()}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"error: {decrypted_message.transpose()}")




