import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpq
import os
import subprocess
import time
import sympy as sp
import ast
from GGH_crypto import GGHCryptosystem, Utils

def embedding_nguyen(basis, point):
        
    if basis.ncols() != point.ncols():
            raise ValueError(f"[Utils] Point is a {1}x{point.ncols()} matrix, but basis is a {basis.nrows()}x{basis.ncols()} one")
    
    n = basis.nrows()

    # Create the initial matrix with R and a column of zeros
    matrix_emb = fmpz_mat([[int(basis[i,j]) if j < n else 0 for j in range(n+1)] for i in range(n)])

    # Add t as the last row, with 1 as the last element
    last_row = [int(point[0,i]) for i in range(n)] + [1]
    matrix_emb = fmpz_mat(matrix_emb.tolist() + [last_row])
    
    matrix_emb, _ = Utils.BKZ_reduction(matrix_emb, block=20, precision=100)

    for i in range(n + 1):
        i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n + 1)]])
        all_elements_are_one = True
        
        for j in range(n + 1):
            if abs(int(i_th_row[0, j])) != 1:
                all_elements_are_one = False
        if all_elements_are_one:
            shortest_vector = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])
            break
    
    if all_elements_are_one:
        print(shortest_vector)
        return shortest_vector

    while True:
        matrix_emb, _ = Utils.BKZ_reduction(matrix_emb, 60, pruned=True, precision=100, bkzautoabort=False, bkzmaxloops=5, nolll=True)

        for i in range(n + 1):
            i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n + 1)]])
            all_elements_are_one = True
            for j in range(n + 1):
                if abs(int(i_th_row[0, j])) != 1:
                    all_elements_are_one = False
            if all_elements_are_one:
                shortest_vector = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])
                break

        if all_elements_are_one:
            print(shortest_vector)
            return shortest_vector

def flint_to_sympy(basis_flint):
    cols = basis_flint.ncols()
    rows = basis_flint.nrows()
    sympy_matrix = sp.zeros(rows, cols)
    for i in range(rows):
        for j in range(cols):
            sympy_matrix[i,j] = sp.Rational(basis_flint[i, j].str())
    return sympy_matrix

def write_matrix_to_file(matrix, filename):
    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)

    rows = matrix.nrows()
    cols = matrix.ncols()

    # Open the file for writing
    with open(filename, "w") as file:
        file.write("[")
        # Iterate over each row of the matrix
        for i in range(rows):
            # Iterate over each column of the matrix
            file.write("[")
            for j in range(cols):
                # Write the element to the file
                file.write(matrix[i, j].str() + " ")
            # Write a newline character after each row
            file.write("]\n")
        file.write("]")
        
    return filename

def matrix_inverse_mod(matrix, n):
    
    det = matrix.det()
    
    det_inverse = pow(det, -1, n)
    
    adjugate = det_inverse * det * matrix.inv()

    inverse = modulo_fmpz_mat(adjugate, n)

    return inverse

def modulo_fmpz_mat(matrix, mod):
    rows, cols = matrix.nrows(), matrix.ncols()
    result = fmpz_mat(rows, cols)
    for i in range(rows):
        for j in range(cols):
            result[i,j] = int(matrix[i,j]) % mod
            
    return result

def Nguyen(public_basis, sigma, ciphertext, message):
    
    B = public_basis
    c = ciphertext
    
    mod = 2 * sigma
        
    cs = fmpq_mat(c.nrows(), c.ncols())
    
    # Creazione della matrice cs aggiungendo sigma agli elementi di ciphertext
    for i in range(c.nrows()):
        for j in range(c.ncols()):
            cs[i, j] = c[i, j] + sigma

    # Calcolo dell'inversa di B modulo 2 * sigma
    B_inv_mod = matrix_inverse_mod(B, mod)

    # Prima stima di m_2sigma
    m_2sigma = cs * B_inv_mod
    m_2sigma = fmpq_mat(modulo_fmpz_mat(m_2sigma, mod))
    
    mod = fmpq(mod)
    
    # Miglioramento del CVP usando la tecnica di Nguyen
    better_CVP = 2 * (fmpq_mat(c - (m_2sigma * B)) / mod)  # 2 * better_CVP, pagina 10 di Nguyen
    
    # Calcolo dell'errore trovato
    error_found = fmpq_mat(embedding_nguyen(B, better_CVP))

    # Calcolo di m_prime
    m_prime = better_CVP - error_found
    B_inv = (2 * B).inv()  # Inversa di 2 * B
    m_prime = m_prime * B_inv

    # Calcolo del risultato finale
    final_result = m_2sigma + (mod * m_prime)

    # Verifica del risultato
    if final_result != message:
        print("wrong sign.. oopsie")
        # Se il risultato non è corretto, inverte il segno di error_found e ricalcola
        error_found = -error_found
        m_prime = better_CVP - error_found
        m_prime = m_prime * B_inv
        final_result = m_2sigma + (mod * m_prime)

    return final_result
    
dimension = 300
tries = 0
print("Finding a basis which can be inverted mod 6...")
while True:
    GGH_object = GGHCryptosystem(dimension = dimension)
    
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

print(f"time taken: {time.time() - start_time}")

print(decrypted_message == message)

print(GGH_object.error)