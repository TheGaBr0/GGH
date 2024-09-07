from GGH_crypto import Utils, GGHCryptosystem
from flint import fmpz_mat
from decimal import Decimal
import time
import numpy as np

def gs(basis):
    ortho_basis = basis[0:1,:].copy()
    for i in range(1, basis.shape[0]):
        proj = np.diag((basis[i,:].dot(ortho_basis.T)/np.linalg.norm(ortho_basis,axis=1)**2).flat).dot(ortho_basis)
        ortho_basis = np.vstack((ortho_basis, basis[i,:] - proj.sum(0)))
    return ortho_basis

def min_norm_row(matrix):
    norms = []
    for j in range(matrix.nrows()):
        row = [matrix[j, i] for i in range(matrix.ncols())]
        norms.append(Utils.vector_l2_norm(row))
    min_norm = min(norms)
    return norms[norms.index(min_norm)]

def calculate_rho(basis, numpy=False):

    if numpy:
        basis_orthogonalized = gs(np.array(basis.tolist()).astype(int))
        basis_orthogonalized = Utils.npsp_to_fmpq_mat(basis_orthogonalized)
    else:
        basis_orthogonalized = Utils.gram_schmidt(basis)
    min_norm = min_norm_row(basis_orthogonalized)
    rho = Decimal(0.5) * min_norm

    print(f'[GGH-HNF] Rho is: {rho:f}')
    return rho

dimension = 50

GGH_object = GGHCryptosystem(dimension = dimension)
GGH_object.encrypt()

time_now = time.time()
print(calculate_rho(GGH_object.private_basis, numpy=False))
print(time.time()-time_now)

time_now = time.time()
print(calculate_rho(GGH_object.private_basis, numpy=True))
print(time.time()-time_now)

