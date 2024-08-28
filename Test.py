from mpmath import (matrix, ones, zeros, randmatrix, nprint, chop, iv,
                    lu_solve, residual, fp, lu, diag, eye, eps, qr, mpf, norm)
from flint import fmpz_mat, fmpz, fmpq, fmpq_mat
import math
import random
import sympy as sp
from GGH_crypto import Utils
import time
from decimal import Decimal, getcontext

dimension = 100000
random_elements = [random.randint(-128, 127) for _ in range(dimension)]
message = fmpq_mat([random_elements])

dimension = 3

time_start = time.time()
print(max(sum(abs(Decimal(int(x.numer())) / Decimal(int(x.denom()))) for x in row) for row in message.tolist()))
print(time.time()-time_start)

m_2 = matrix([[item.str() for item in sublist] for sublist in message.tolist()])
time_start = time.time()
print(norm(m_2, p=1))
print(time.time()-time_start)