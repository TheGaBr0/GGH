from flint import *
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from decimal import Decimal
import math
import random
from LLL import LLL_reduction

def embedding(R, t):
    n = len(t)
    
    matrix = np.zeros((n+1, n+1))

    matrix[0:n, 0] = t

    matrix[0:n, 1:n+1] = R
    
    matrix[n, 0] = 1
    
    matrix = LLL_reduction(matrix.T, n+1).T #my implementation of LLL is based on row vectors matrices

    first_column = matrix[:-1, 0]
    
    
    return t - first_column

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
    
    plt.scatter(t_np[0],t_np[1])
    
    plt.scatter(int(w_2[0]),int(w_2[1]), color='red', s=70)
    
    plt.scatter(int(w[0]),int(w[1]), color='blueviolet', s=70)
    
    plt.annotate(label_t, # this is the text
                 (t_np[0],t_np[1]), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,-12), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
    plt.annotate(label_w, 
                 (int(w[0]),int(w[1])), 
                 textcoords="offset points",
                 xytext=(0,10), 
                 ha='center') 
    plt.annotate(label_w_2, 
                 (int(w_2[0]),int(w_2[1])), 
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
    
def column_norm(col):
    return Decimal(sum(Decimal(fmpq_to_decimal(x))**2 for x in col)).sqrt()
    
def get_hadamard_ratio(basis):
    norms = []
        
    for j in range(basis.ncols()):
        column = [basis[i, j] for i in range(basis.nrows())]
        norms.append(Decimal(sum(Decimal(int(x))**2 for x in column)).sqrt())
    
    denominator = math.prod(norms)
    numerator = abs(Decimal(basis.det().str()))
    result = (numerator / denominator) ** Decimal(1 / 2)
    return f"{result:.16f}"


  
R = fmpz_mat([[1,2],[3,0]])


B = fmpz_mat([[12,-6],[7,-4]])

print(f"Good basis R = {R.tolist()} with det = {R.det()} and Hadamard ratio = {get_hadamard_ratio(R)}")
print(f"Bad basis B = {B.tolist()} with det = {B.det()} and Hadamard ratio = {get_hadamard_ratio(B)}")

R_np = np.array(R.tolist()).astype(int)
B_np = np.array(B.tolist()).astype(int)
t_np = np.array([6, 3]).astype(int)

w = embedding(R_np.T, t_np.T)
w_2 = embedding(B_np.T, t_np.T)
print(f"CVP found by R is {w.tolist()}, t-w = {column_norm(numpy_to_fmpz_mat([t_np-w]))}")
print(f"CVP found by B is {w_2.tolist()}, t-w_2 = {column_norm(numpy_to_fmpz_mat([t_np-w]))}")

generate_lattice_points(R_np, B_np, w, w_2, t_np)




