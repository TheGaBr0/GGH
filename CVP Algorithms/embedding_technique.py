from flint import *
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from decimal import Decimal
import math
import random
from GGH import GGHCryptosystem

def embedding(R, t):
    n = R.nrows()

    # Create the initial matrix with R and a column of zeros
    matrix_emb = fmpz_mat([[int(R[i,j]) if j < n else 0 for j in range(n+1)] for i in range(n)])

    # Add t as the last row, with 1 as the last element
    last_row = [int(t[0,i]) for i in range(n)] + [1]
    matrix_emb = fmpz_mat(matrix_emb.tolist() + [last_row])
    
    matrix_emb = matrix_emb.lll()

    print(matrix_emb)

    min_norm = np.inf
    shortest_row = t

    for i in range(n+1):
        i_th_row = fmpz_mat([[matrix_emb[i, j] for j in range(n+1)]])
        norm = row_norm(i_th_row)
        if norm < min_norm:
            min_norm = norm
            # We take the first n values of the shortest row
            shortest_row = fmpz_mat([[matrix_emb[i, j] for j in range(n)]])
    
    return t - shortest_row

def numpy_to_fmpz_mat(numpy_matrix):
        return fmpz_mat([[int(item) for item in sublist] for sublist in numpy_matrix])

def generate_lattice_points(R_np, B_np, w, w_2, t_np, limit=5):
    # Create a meshgrid of integer coordinates
    x = np.arange(-limit, limit + 1)
    y = np.arange(-limit, limit + 1)
    coords = np.array(np.meshgrid(x, y)).T.reshape(-1, 2)

    # Multiply the meshgrid by the basis matrix to get lattice points
    lattice_points = np.dot(coords, R_np)
    
    plt.plot(lattice_points[:, 0], lattice_points[:, 1], 'bo', color='black')
    plt.grid(True)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Lattice Plot')
    label_w = "w"
    label_w_2 = "w_2"
    label_t = "t"
    
    plt.scatter(t_np[0][0],t_np[0][1])
    
    plt.scatter(int(w_2[0, 0]),int(w_2[0, 1]), color='red', s=70)
    
    plt.scatter(int(w[0, 0]),int(w[0, 1]), color='blueviolet', s=70)
    
    plt.annotate(label_t, # this is the text
                 (t_np[0][0],t_np[0][1]), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,-12), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
    plt.annotate(label_w, 
                 (int(w[0, 0]),int(w[0, 1])), 
                 textcoords="offset points",
                 xytext=(0,10), 
                 ha='center') 
    plt.annotate(label_w_2, 
                 (int(w_2[0, 0]),int(w_2[0, 1])), 
                 textcoords="offset points", 
                 xytext=(0,-12), 
                 ha='center') 
    
    
    plt.arrow(0, 0, R_np[0][0], R_np[0][1], lw=2, head_width=0.1, head_length=0.1, fc='blueviolet', ec='blueviolet')
    plt.arrow(0, 0, R_np[1][0], R_np[1][1], lw=2, head_width=0.1, head_length=0.1, fc='blueviolet', ec='blueviolet')
    plt.arrow(0, 0, B_np[0][0], B_np[0][1], lw=2, head_width=0.1, head_length=0.1, fc='red', ec='red')
    plt.arrow(0, 0, B_np[1][0], B_np[1][1], lw=2, head_width=0.1, head_length=0.1, fc='red', ec='red')  
    
    plt.show()
    
def fmpq_to_decimal(fmpq_number):
    fraction_str = fmpq_number.str()
    try:
        if '/' in fraction_str:
            numerator_str, denominator_str = fraction_str.split('/')
            numerator = Decimal(numerator_str)
            denominator = Decimal(denominator_str)
            result = numerator / denominator
        else:
            result = Decimal(fraction_str)
        return result
    except Exception as e:
        print(f"Error converting {fraction_str} to decimal: {e}")
        return None
    
def row_norm(row):
    return Decimal(sum(Decimal(fmpq_to_decimal(x))**2 for x in row)).sqrt()
    
def get_hadamard_ratio(basis):
    norms = []
        
    for i in range(basis.nrows()):
        row = [basis[i, j] for j in range(basis.ncols())]
        norms.append(Decimal(sum(Decimal(int(x))**2 for x in row)).sqrt())
    
    denominator = math.prod(norms)
    numerator = abs(Decimal(basis.det().str()))
    result = (numerator / denominator) ** Decimal(1 / 2)
    return f"{result:.16f}"


GGH_object = GGHCryptosystem(dimension = 2)
GGH_object.encrypt()

R = fmpz_mat([[1,2],[3,0]])

B = fmpz_mat([[5, 4], [-6, -6]])
T = fmpq_mat([[5, 3]])


print(f"Good basis R = {R.tolist()} with det = {R.det()} and Hadamard ratio = {get_hadamard_ratio(R)}")
print(f"Bad basis B = {B.tolist()} with det = {B.det()} and Hadamard ratio = {get_hadamard_ratio(B)}")

R_np = np.array(R.tolist()).astype(int)
B_np = np.array(B.tolist()).astype(int)
t_np = np.array(T.tolist()).astype(int)

w = embedding(R, T)
w_np = np.array(w.tolist()).astype(int)
w_2 = embedding(B, T)
w_2_np = np.array(w_2.tolist()).astype(int)

print(f"CVP found by R is {w.tolist()}, t-w = {row_norm(numpy_to_fmpz_mat([T-w]))}")
print(f"CVP found by B is {w_2.tolist()}, t-w_2 = {row_norm(numpy_to_fmpz_mat([T-w]))}")

generate_lattice_points(R_np, B_np, w_np, w_2_np, t_np)