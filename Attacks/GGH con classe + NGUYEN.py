from GGH import GGHCryptosystem
import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpz
import os
import subprocess
import time
import sympy as sp
import ast

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

def embedding(B, c):
    # Aggiungi una colonna di zeri a destra di B

    n = B.ncols()

    matrix_emb = fmpz_mat([[B[i,j] if i < n and j < n else 0 for j in range(n+1)] for i in range(n+1)])
    
    t_to_fmpz = fmpz_mat([[int(x) for x in c]])

    for j in range(n):
        matrix_emb[n, j] = t_to_fmpz[0, j]

    matrix_emb[n,n] = 1

    return matrix_emb

def reduction(matrix, block = 20, pruned=False):
     
    path = write_matrix_to_file(matrix, 'BE.txt')

    start_time2 = time.time()
    if os.name == 'nt':
        
        path = '/mnt/' + path.replace(":", "") # removing ':' from windows path (c:\Users...) because wsl uses this format: mnt/c/Users...
        path = path.replace("\\", "/") # reversing \ with /
        path = path.replace(" ", "\ ") # fixing folders with spaces in the name
        
        if not pruned:
            command = f"wsl fplll {path} -a bkz -b {block} -p 85 -f mpfr -bkzautoabort > out.txt" 
        else:
            command = f"wsl fplll {path} -a bkz -b {block} -p 85 -f mpfr -s default.json -bkzmaxloops 3 > out.txt" 
    else:
        command = f"fplll '{path}' -a bkz -b 20 -p 100 -f mpfr -m proved "

    try:
        # Run the command and capture its output
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = process.communicate()
        output = output.decode('utf-8')
        error = error.decode('utf-8')

        # Check if there was an error
        if error:
            print("Error:", error)
       
        print(f"reduction completed, time taken: {time.time() - start_time2}")

        return load_matrix_from_file("out.txt", True)

    except Exception as e:
        return None, str(e)
    
def load_matrix_from_file(filename, fplll = False):
    
    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)

    if not fplll:

        with open(filename, 'r') as file:
            matrix_data = file.read().strip()

        # Removing the outer brackets
        matrix_data = matrix_data[1:-1]

        # Splitting into rows and then into individual elements
        rows = matrix_data.split('\n')
        matrix_list = [list(map(int, row.strip()[1:-1].split())) for row in rows]
        print(matrix_list)
        return fmpz_mat(matrix_list)
    else:
        with open(filename, 'r') as file:
            content = file.read()
        
        # Remove whitespace and newlines
        content = content.replace(']', '],')
        content = content.replace(' ]', ']')
        content = content.replace(' ', ', ')
        content = content[:-5] + ']'

        # Use ast.literal_eval to safely evaluate the string as a Python expression
        matrix = ast.literal_eval(content)
        
        return fmpz_mat(matrix)


def modulo_fmpz_mat(matrix, mod):
    rows, cols = matrix.nrows(), matrix.ncols()
    result = fmpz_mat(rows, cols)
    for i in range(rows):
        for j in range(cols):
            result[i,j] = int(matrix[i,j]) % mod
            
    return result

def get_first_row_no_last_element(matrix):
    cols = matrix.ncols()  # Numero di colonne
    first_row_no_last_element = [matrix[0, j] for j in range(cols-1)]
    return fmpz_mat([first_row_no_last_element])

def Nguyen(public_basis, sigma, ciphertext):
    
    B = public_basis
    c = ciphertext
    
    mod = 2*sigma
    
    cs = fmpq_mat(c.nrows(), c.ncols())
    
    for i in range(c.nrows()):
        for j in range(c.ncols()):
            cs[i,j] = c[i,j] + sigma

    B_inv_mod = matrix_inverse_mod(B, mod)

    m_2sigma = cs * B_inv_mod
    
    m_2sigma = modulo_fmpz_mat(m_2sigma, mod)

    write_matrix_to_file(m_2sigma, "m_2sigma.txt")
    
    better_CVP = 2 * (fmpq_mat(c - (m_2sigma * B)) / mod) #2 * better_CVP, pagina 10 nguyen, in fondo

    
    write_matrix_to_file(better_CVP, "better_CVP.txt")

    BE = embedding(B, better_CVP)

    print("embedding step completed")

    reduced_BKZ = reduction(reduction(BE), 60, True)

    write_matrix_to_file(reduced_BKZ, "reduced2.txt")

    print("reduction step completed")

    e = get_first_row_no_last_element(reduced_BKZ)

    write_matrix_to_file(e, "e.txt")

    B_inv = (2*B).inv()

    B_inv = flint_to_sympy(B_inv)
    
    better_CVP = flint_to_sympy(better_CVP)

    e = flint_to_sympy(e)

    m_2sigma = flint_to_sympy(m_2sigma)

    m_prime = better_CVP - e

    m_prime = m_prime * B_inv

    print(m_prime)

    final_result = m_2sigma + (mod*m_prime)

    return final_result
    
dimension = 3

B = fmpz_mat([[145, -73, -23],[-39,  21,  16],[-165,  80,  11]])
c = fmpz_mat([[4452, -1964, 735]])
sigma = 3

#print(f"message: {message.transpose()}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"sigma: {sigma}")

start_time = time.time()

decrypted_message = Nguyen(B, sigma, c)

print(f"time taken: {time.time() - start_time}")

#print(f"decrypted message: {decrypted_message}")

print(decrypted_message)





