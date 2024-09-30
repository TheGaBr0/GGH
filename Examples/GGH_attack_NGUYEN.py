from flint import fmpz_mat, fmpq_mat, fmpq
import time
import sympy as sp
from GGH_crypto import GGHCryptosystem

def embedding_nguyen(basis, point):
    """
    Implement an iterative embedding technique for lattice-based cryptanalysis. Due to the 
    limitations of FPLLL, the attack is executed iteratively, with increasing reduction strength in each iteration.
    
    :param basis: The lattice basis
    :param point: The target point
    :return: The shortest vector found in form of an array with +-1 as elements
    """
    if basis.ncols() != point.ncols():
        raise ValueError(f"[Utils] Point is a {1}x{point.ncols()} matrix, but basis is a {basis.nrows()}x{basis.ncols()} one")
    
    n = basis.nrows()

    # Create the initial matrix with R and a column of zeros
    matrix_emb = fmpz_mat([[int(basis[i,j]) if j < n else 0 for j in range(n+1)] for i in range(n)])

    # Add t as the last row, with 1 as the last element
    last_row = [int(point[0,i]) for i in range(n)] + [1]
    matrix_emb = fmpz_mat(matrix_emb.tolist() + [last_row])
    
    # Perform initial BKZ reduction
    matrix_emb, _ = Utils.BKZ_reduction(matrix_emb, block=20, precision=100)

    # Check for a vector with all elements ±1
    for i in range(n + 1):
        i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n + 1)]])
        if all(abs(int(i_th_row[0, j])) == 1 for j in range(n + 1)):
            return fmpz_mat([[matrix_emb[i, j] for j in range(n)]])

    # If no suitable vector found, continue with stronger BKZ reduction
    # Each 5 loops we'll check if BKZ has found the vector we're looking for 
    while True:
        matrix_emb, _ = Utils.BKZ_reduction(matrix_emb, 60, pruned=True, precision=100, bkzautoabort=False, bkzmaxloops=5, nolll=True)

        for i in range(n + 1):
            i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n + 1)]])
            if all(abs(int(i_th_row[0, j])) == 1 for j in range(n + 1)):
                return fmpz_mat([[matrix_emb[i, j] for j in range(n)]])

def matrix_inverse_mod(matrix, n):
    """
    Compute the modular inverse of a matrix.
    
    :param matrix: The input matrix
    :param n: The modulus
    :return: The modular inverse of the matrix
    """
    det = matrix.det()
    det_inverse = pow(det, -1, n)
    adjugate = det_inverse * det * matrix.inv()
    return modulo_fmpz_mat(adjugate, n)

def modulo_fmpz_mat(matrix, mod):
    """
    Apply modulo operation to each element of the matrix.
    
    :param matrix: The input matrix
    :param mod: The modulus
    :return: The resulting matrix after modulo operation
    """
    rows, cols = matrix.nrows(), matrix.ncols()
    result = fmpz_mat(rows, cols)
    for i in range(rows):
        for j in range(cols):
            result[i,j] = int(matrix[i,j]) % mod
    return result

def Nguyen(public_basis, sigma, ciphertext, message):
    """
    Implement Nguyen's attack on the GGH cryptosystem.
    
    :param public_basis: The public basis of the GGH cryptosystem
    :param sigma: The sigma parameter of the GGH cryptosystem
    :param ciphertext: The ciphertext to be decrypted
    :param message: The original message (for verification)
    :return: The decrypted message
    """
    B = public_basis
    c = ciphertext
    mod = 2 * sigma
    
    # Create cs by adding sigma to each element of the ciphertext
    cs = fmpq_mat(c.nrows(), c.ncols())
    for i in range(c.nrows()):
        for j in range(c.ncols()):
            cs[i, j] = c[i, j] + sigma

    # Calculate the inverse of B modulo 2 * sigma
    B_inv_mod = matrix_inverse_mod(B, mod)

    # First estimate of m_2sigma
    m_2sigma = cs * B_inv_mod
    m_2sigma = fmpq_mat(modulo_fmpz_mat(m_2sigma, mod))
    
    #everything must be converted to fmpq to avoid FLINT arithmetic issues. This is just to be safe
    mod = fmpq(mod)
    
    # Improve CVP using Nguyen's technique
    better_CVP = 2 * (fmpq_mat(c - (m_2sigma * B)) / mod)
    
    # Calculate the found error
    error_found = fmpq_mat(embedding_nguyen(B, better_CVP))

    # Calculate m_prime
    m_prime = better_CVP - error_found
    B_inv = (2 * B).inv()
    m_prime = m_prime * B_inv

    # Calculate the final result
    final_result = m_2sigma + (mod * m_prime)

    # Verify the result
    if final_result != message:
        print("Solution computed is wrong, let's try to change its sign...")
        # If the result is incorrect, invert the sign of error_found and recalculate
        error_found = -error_found
        m_prime = better_CVP - error_found
        m_prime = m_prime * B_inv
        final_result = m_2sigma + (mod * m_prime)

    return final_result

dimension = 100
tries = 0
print("Finding a basis which can be inverted mod 6...")

#The simplest case is when the public basis is invertible mod 2*sigma
#Therefore, we will generate GGH instances repeatedly until we find a matrix that satisfies this invertibility condition
while True:
    GGH_object = GGHCryptosystem(dimension=dimension)
    
    B, sigma = GGH_object.public_key
    if sp.gcd(int(B.det()), 2*sigma) == 1:
        GGH_object.encrypt()
        print(f"Try number {tries}: success")
        break
    else:
        tries += 1
        print(f"Try number {tries}: fail")

message = GGH_object.message
print(f"message: {message}")

start_time = time.time()
decrypted_message = Nguyen(B, sigma, GGH_object.ciphertext, message)
print(f"Total time taken to attack: {time.time() - start_time}")

if decrypted_message == message:
    print(f"Attack completed with success, decrypted message:\n{decrypted_message}")