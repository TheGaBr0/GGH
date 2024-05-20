from GGH import GGHCryptosystem
import sympy as sp
import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpz
import os
import subprocess
import time
dimension = 2
GGH_object = GGHCryptosystem(dimension = dimension)
GGH_object.encrypt()

R = GGH_object.private_key[1]
B = GGH_object.public_key[0]

print(f"R: {R}")

print(f"B: {B}")
    
message = GGH_object.message

ciphertext = GGH_object.ciphertext

print(f"C: {ciphertext}")

decrypted_message = GGH_object.decrypt()

print(f"sigma: {GGH_object.public_key[1]}")
print(f"message: {message.transpose()}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"decrypted_message: {decrypted_message.transpose()}")




