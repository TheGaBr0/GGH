from GGH_crypto import GGHHNFCryptosystem, Utils
from flint import fmpz_mat
import time

# Set the dimension for the cryptosystem
dimension = 100

# Initialize the GGH cryptosystem
GGH_object = GGHHNFCryptosystem(dimension=dimension, GGH_private=False, alpha=0.75, debug=True)

# Encrypt a message
GGH_object.encrypt()

# Extract the public key, private key, error vector, and ciphertext
B = GGH_object.public_key[0]
R = GGH_object.private_key
error = GGH_object.error
ciphertext = GGH_object.ciphertext

decrypted_message = GGH_object.decrypt()

# Print the error vector and decrypted message

if error != decrypted_message:
    print("Decryption failed, try with GGH_private = True or lower the alpha value")
    exit()
print(f"error: {error}")
print(f"error: {decrypted_message}")

# Start of the attack using embedding technique
n = B.nrows()

# Create the initial matrix with B (public key) and a column of zeros
matrix_emb = fmpz_mat([[int(B[i,j]) if j < n else 0 for j in range(n+1)] for i in range(n)])

# Add ciphertext as the last row, with 1 as the last element
last_row = [int(ciphertext[0,i]) for i in range(n)] + [1]
matrix_emb = fmpz_mat(matrix_emb.tolist() + [last_row])

# Start timing the attack
time_now = time.time()
found = False

# First attempt with lower block size
while True:
    # Perform BKZ reduction on the matrix
    matrix_emb, error_message = Utils.BKZ_reduction(matrix_emb, block=20, pruned=True, precision=100, bkzmaxloops=5)

    # Extract the first row of the reduced matrix
    i_th_row = fmpz_mat([[matrix_emb[0, j] for j in range(n)]])

    # Check if the extracted row matches the error vector
    if error == i_th_row or error == -1*i_th_row or error_message is None:
        if error == i_th_row or error == -1*i_th_row:
            found = True
            print(f"error: {i_th_row}")
        break
    else:
        print(error_message)

# If not found, try again with higher block size
if not found:
    while True:
        # Perform BKZ reduction with increased block size
        matrix_emb, error_message = Utils.BKZ_reduction(matrix_emb, precision=100, block=60, pruned=True, nolll=True, bkzmaxloops=5)

        # Extract the first row of the reduced matrix
        i_th_row = fmpz_mat([[matrix_emb[0, j] for j in range(n)]])

        # Check if the extracted row matches the error vector
        if error == i_th_row or error == -1*i_th_row or error_message is None:
            print(f"error: {i_th_row}")
            break

print(f"Total time taken to attack: {time.time() - time_now}")
