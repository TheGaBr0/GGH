import numpy as np
import random
from copy import deepcopy
import math

def gram_schmidt(matrix):
    n = matrix.shape[0]
    
    w_array = [deepcopy(matrix[0].astype(float))]
    
    for i in range(1, n):
        v_i = matrix[i].astype(float)
        w_i = v_i
        
        for j in range(i):
            w_j = deepcopy(w_array[j])
            
            w_i -= np.dot(v_i, w_j) / np.dot(w_j, w_j) * w_j
        w_array += [w_i]
    
    #When you are using lattice basis reduction algorithms like the 
    # LLL (Lenstra-Lenstra-Lovász) algorithm to reduce a lattice basis, 
    # you do not need to explicitly normalize the orthogonalized vectors 
    # as you would in the standard Gram-Schmidt process. 
    # The LLL algorithm inherently takes care of the orthogonalization and size reduction of the basis vectors.

    return w_array


def mu(b_i, b_j):
    return np.dot(b_i, b_j) / np.dot(b_j, b_j) if np.dot(b_j, b_j) != 0 else 0

def LLL_reduction(b, n):
    
    k = 1
    
    d = 0.75
    
    basis = deepcopy(b)
    
    b_prime = np.array(gram_schmidt(basis))
    
    while k < n:
        for j in range(k - 1, -1, -1): #da k-1 a 0 con step -1
            mu_value = mu(basis[k], b_prime[j])
            if abs(mu_value) > 0.5:
                basis[k] = basis[k] - round(mu_value) * basis[j]
                b_prime = np.array(gram_schmidt(basis))

        
        if np.dot(b_prime[k], b_prime[k]) >= (d - mu(basis[k], b_prime[k-1])**2)* (np.dot(b_prime[k-1], b_prime[k-1])):
            k += 1
        else:
            basis[[k, k-1]] = basis[[k-1, k]]
            b_prime = np.array(gram_schmidt(basis))
            k = max(k - 1, 1)
           
    return basis


