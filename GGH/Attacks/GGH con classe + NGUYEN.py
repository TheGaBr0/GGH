from GGH import GGHCryptosystem
import sympy as sp
import numpy as np
from flint import fmpz_mat, fmpq_mat, fmpz
import os
import subprocess
import symengine as sym
import time

def write_matrix_to_file(matrix, filename):
    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
    
    if isinstance(matrix, sp.matrices.dense.MutableDenseMatrix):
        rows = matrix.shape[0]
        cols = matrix.shape[1]
    else:
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
                file.write(str(matrix[i, j]) + " ")
            # Write a newline character after each row
            file.write("]\n")
        file.write("]")
        
    return filename

def sympy_to_rational_flint(basis_sympy):
    rows, cols = basis_sympy.shape
    flint_matrix = fmpq_mat(rows, cols)
    for i in range(rows):
        for j in range(cols):
            fraction_val = sp.fraction(basis_sympy[i, j])
            flint_matrix[i,j] = fmpq(int(fraction_val[0]), int(fraction_val[1]))
    return flint_matrix

def sympy_to_fmpz_mat(basis_sympy):
    return fmpz_mat([[int(item) for item in sublist] for sublist in basis_sympy.tolist()])

def matrix_inverse_mod(matrix, n):
    
    det = matrix.det()
    
    det_inverse = pow(det, -1, n)
    
    adjugate = det_inverse * det * matrix.inv()

    inverse = modulo_fmpz_mat(adjugate, n)

    return inverse

def embedding(B, c):
    # Aggiungi una colonna di zeri a destra di B
        
    B_sympy = flint_to_sympy(B)
    
    c_sympy = flint_to_sympy(c)
        
    B_sympy = B_sympy.row_join(sp.Matrix.zeros(B_sympy.rows, 1))
    
    # Aggiungi una riga con il vettore c e un 1 come ultimo elemento
    c_sympy = c_sympy.col_insert(c_sympy.cols, sp.Matrix([1]))
    B_sympy = B_sympy.row_insert(B_sympy.rows, c_sympy)
    
    result = sympy_to_fmpz_mat(B_sympy)
    return result

def reduction(matrix):        
    path = write_matrix_to_file(matrix, 'BE.txt')
    
    if os.name == 'nt':
        
        path = '/mnt/' + path.replace(":", "") # removing ':' from windows path (c:\Users...) because wsl uses this format: mnt/c/Users...
        path = path.replace("\\", "/") # reversing \ with /
        path = path.replace(" ", "\ ") # fixing folders with spaces in the name
        
        command = f"wsl fplll {path} -a bkz -b 20 -f mpfr -p 70"
    else:
        command = f"fplll '{path}' -a bkz -b 20 -f mpfr -p 70"

    try:
        # Run the command and capture its output
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, error = process.communicate()
        output = output.decode('utf-8')
        error = error.decode('utf-8')

        # Check if there was an error
        if error:
            print("Error:", error)
        else:
            output = output.replace(']', '],')
            output = output.replace(' ]', ']')
            output = output.replace(' ', ', ')
            output = output[:-5] + ']'

            return fmpz_mat(eval(output))

    except Exception as e:
        return None, str(e)
    
def load_matrix_from_file(filename):
    
    filename = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)

    with open(filename, 'r') as file:
        matrix_data = file.read().strip()

    # Removing the outer brackets
    matrix_data = matrix_data[1:-1]

    # Splitting into rows and then into individual elements
    rows = matrix_data.split('\n')
    matrix_list = [list(map(int, row.strip()[1:-1].split())) for row in rows]

    # Creating a SymPy matrix
    return sympy_to_fmpz_mat(sp.Matrix(matrix_list))

def flint_to_sympy(basis_flint):
    cols = basis_flint.ncols()
    rows = basis_flint.nrows()
    sympy_matrix = sp.zeros(rows, cols)
    for i in range(rows):
        for j in range(cols):
            sympy_matrix[i,j] = sp.Rational(basis_flint[i, j].str())
    return sympy_matrix


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
    
    B = public_basis.transpose()
    c = ciphertext.transpose()
    
    mod = 2*sigma
    
    cs = fmpq_mat(c.nrows(), c.ncols())
    
    for i in range(c.nrows()):
        for j in range(c.ncols()):
            cs[i,j] = c[i,j] + sigma
            
    
    B_inv_mod = matrix_inverse_mod(B, mod)

    m_2sigma = cs * B_inv_mod
    
    m_2sigma = modulo_fmpz_mat(m_2sigma, mod)
    
    better_CVP = 2 * (fmpq_mat(c - (m_2sigma * B)) / mod) #2 * better_CVP, pagina 10 nguyen, in fondo
    
    BE = embedding(B, better_CVP)
    
    # Run the command
    reduced_BKZ = reduction(BE)

    #reduced_BKZ = optimized_LLL(BE)

    e = get_first_row_no_last_element(reduced_BKZ)
    
    m_prime = (better_CVP - e) * (2 * B).inv()
    
    write_matrix_to_file(m_prime, "m_prime.txt")
    
    B = flint_to_sympy(B)
    
    B = 2 * B
    
    better_CVP = flint_to_sympy(better_CVP)
    
    e = flint_to_sympy(e)
    
    B_inv = B._rep.to_field().inv().to_Matrix()
    
    m_prime = (better_CVP - e) * B_inv
    
    write_matrix_to_file(m_prime, "m_prime2.txt")
    
    m_2sigma = flint_to_sympy(m_2sigma)
    
    final_result = m_2sigma + (mod*m_prime)

    return final_result
    
dimension = 50
tries = 0
while True:
    GGH_object = GGHCryptosystem(dimension = dimension)
    GGH_object.encrypt()
    
    B, sigma = GGH_object.public_key
    
    if sp.gcd(int(B.det()), 2*sigma) == 1:
        break
    else:
        tries += 1
        print(tries)
    
message = GGH_object.message

ciphertext = GGH_object.ciphertext

decrypted_message = GGH_object.decrypt()

print(f"message: {message.transpose()}")
# print(f"decrypted message: {decrypted_message.transpose()}")
print(f"sigma: {sigma}")

start_time = time.time()

decrypted_message = Nguyen(B, sigma, ciphertext)

print(f"time taken: {time.time() - start_time}")

print(f"decrypted message: {decrypted_message}")

print(sympy_to_fmpz_mat(decrypted_message).transpose() == message)




