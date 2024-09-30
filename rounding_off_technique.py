from flint import fmpq_mat
from GGH_crypto import GGHCryptosystem, Utils

# Initialize the GGH cryptosystem with a 2-dimensional lattice
GGH_object = GGHCryptosystem(dimension=2)

# Extract the private (good) and public (potentially bad) bases
R = GGH_object.private_basis 
B = GGH_object.public_basis   

# Define a target point, possibly not in the lattice
T = fmpq_mat([[12, 21]])

# Print information about both bases
print(f"Good basis R = {R.tolist()} with det = {R.det()} and Hadamard ratio = {Utils.get_hadamard_ratio(R)[1]}")
print(f"Bad basis B = {B.tolist()} with det = {B.det()} and Hadamard ratio = {Utils.get_hadamard_ratio(B)[1]}")

# Solve the Closest Vector Problem (CVP) using Babai's rounding algorithm with the private basis
# This should give a reasonably accurate result
w = Utils.babai_rounding(R, T, visualize=True)

# Solve CVP using Babai's rounding algorithm with the public basis
# This is more likely to give an incorrect result compared to using the private basis
w_2 = Utils.babai_rounding(B, T, visualize=True)

# Print the results and the distances between the target point and the found lattice points
print(f"CVP found by R is {w.tolist()}, t-w = {Utils.vector_l2_norm(T-w)}")
print(f"CVP found by B is {w_2.tolist()}, t-w_2 = {Utils.vector_l2_norm(T-w_2)}")

# Note: Babai's rounding algorithm is an approximate method for solving CVP.
# It generally performs well with a good basis like R, but may produce
# significant errors with a "bad" basis like B.
# The quality of the result using B heavily depends on how "bad" B is compared to R.
# In particular:
# 1. B is likely to produce more errors than R in most cases.
# 2. The worse the Hadamard ratio of B, the more likely it is to produce incorrect results.
