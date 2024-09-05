from GGH_crypto import Utils, GGHCryptosystem
from flint import fmpz_mat

dimension = 100



GGH_object = GGHCryptosystem(dimension = dimension)
GGH_object.encrypt()

B = GGH_object.public_key

R = GGH_object.private_key[1]
    
message = GGH_object.message

ciphertext = GGH_object.ciphertext

decrypted_message = GGH_object.decrypt()

print(type(decrypted_message))

print(f"error: {message}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"error: {decrypted_message}")

