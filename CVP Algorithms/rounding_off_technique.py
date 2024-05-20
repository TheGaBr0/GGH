from flint import *
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from decimal import Decimal
import math
import random

def rounding(R, w):
    
    
    a = R.inv().transpose() * w
    cols = a.ncols()
    rows = a.nrows()

    for i in range(rows):
        for j in range(cols):
            a[i,j] = round(a[i,j])
            
    result = R.transpose() * a
    print(result)
    return result

def generate_lattice_points(R_np, B_np, w, w_2, limit=5):
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
    
    plt.scatter(X_rational,Y_rational)
    
    plt.scatter(int(w_2[0, 0]),int(w_2[1, 0]), color='red', s=70)
    
    plt.scatter(int(w[0, 0]),int(w[1, 0]), color='blueviolet', s=70)
    
    plt.annotate(label_t, # this is the text
                 (X_rational,Y_rational), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(0,-12), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
    plt.annotate(label_w, # this is the text
                 (int(w[0, 0]),int(w[1, 0])), 
                 textcoords="offset points",
                 xytext=(0,10), 
                 ha='center') 
    plt.annotate(label_w_2, 
                 (int(w_2[0, 0]),int(w_2[1, 0])), 
                 textcoords="offset points", 
                 xytext=(0,-12), 
                 ha='center') 
    
    
    plt.arrow(0, 0, R_np[0][0], R_np[0][1], lw=2, head_width=0.1, head_length=0.1, fc='blueviolet', ec='blueviolet')
    plt.arrow(0, 0, R_np[1][0], R_np[1][1], lw=2, head_width=0.1, head_length=0.1, fc='blueviolet', ec='blueviolet')
    plt.arrow(0, 0, B_np[0][0], B_np[0][1], lw=2, head_width=0.1, head_length=0.1, fc='red', ec='red')
    plt.arrow(0, 0, B_np[1][0], B_np[1][1], lw=2, head_width=0.1, head_length=0.1, fc='red', ec='red')  
    
    plt.show()
    
def get_hadamard_ratio(basis):
    norms = []
        
    for j in range(basis.ncols()):
        column = [basis[i, j] for i in range(basis.nrows())]
        norms.append(Decimal(sum(Decimal(int(x))**2 for x in column)).sqrt())
    
    denominator = math.prod(norms)
    numerator = abs(Decimal(basis.det().str()))
    result = (numerator / denominator) ** Decimal(1 / 2)
    return f"{result:.16f}"

def sympy_to_fmpz_mat(basis_sympy):
        return fmpz_mat([[int(item) for item in sublist] for sublist in basis_sympy.tolist()])
  
R = fmpz_mat([[1,2],[3,0]])


B = fmpz_mat([[12,-6],[7,-4]])

print(f"Good basis B' = {R.tolist()} with det = {R.det()} and Hadamard ratio = {get_hadamard_ratio(R)}")
print(f"Good basis B = {B.tolist()} with det = {B.det()} and Hadamard ratio = {get_hadamard_ratio(B)}")

X_rational = sp.Rational(7)
Y_rational = sp.Rational(3.50)
t = fmpq_mat([[fmpq(X_rational.numerator,X_rational.denominator), fmpq(Y_rational.numerator,Y_rational.denominator)]]).transpose()

R_np = np.array(R.tolist()).astype(int)
B_np = np.array(B.tolist()).astype(int)
w = rounding(R, t)
w_2 = rounding(B, t)

print(f"CVP found by B is {w.tolist()}")
print(f"CVP found by B is {w_2.tolist()}")

generate_lattice_points(R_np, B_np, w, w_2)

