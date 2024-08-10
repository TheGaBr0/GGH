from GGH_HNF import GGHHNFCryptosystem
import sympy as sp
import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpz
import os
import subprocess
import symengine as sym
import time
dimension = 3

R = fmpz_mat([[12, -4, -1], [1,  8, -1], [-4,  1, 14]])

GGH_object = GGHHNFCryptosystem(dimension = dimension, R=R)
GGH_object.encrypt()

H = GGH_object.public_key

R = GGH_object.private_key[1]

print(R)

print(H)
    
print(f"rho R: {GGH_object.calculate_rho(R)}")

print(f"rho H: {GGH_object.calculate_rho(H)}")

print(f"lunghezza e: {GGH_object.vector_norm(GGH_object.error)}")

H = GGH_object.public_key
    
error = GGH_object.error

ciphertext = GGH_object.ciphertext

print(ciphertext)

decrypted_message = GGH_object.decrypt()

print(f"error: {error}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"error: {decrypted_message}")




