from flint import fmpq_mat
from GGH_crypto import GGHCryptosystem, Utils

# Initialize the GGH cryptosystem with a 2-dimensional lattice
GGH_object = GGHCryptosystem(dimension=2)

# Extract the private (good) and public (potentially bad) bases
R = GGH_object.private_basis 
B = GGH_object.public_basis  

# Generate a random point, possibly not in the lattice, with small entries for visualization
T = fmpq_mat([[12, 21]])

# Print information about both bases
print(f"Good basis R = {R.tolist()} with det = {R.det()} and Hadamard ratio = {Utils.get_hadamard_ratio(R)[1]}")
print(f"Bad basis B = {B.tolist()} with det = {B.det()} and Hadamard ratio = {Utils.get_hadamard_ratio(B)[1]}")

# Solve the Closest Vector Problem (CVP) using the embedding technique with the private basis
# This should give an accurate result
w = Utils.embedding(R, T, visualize=True, GGH=True, BKZ=False)

# Solve CVP using the embedding technique with the public basis
# This may give an incorrect result, especially if the Hadamard ratio of B is low (<0.50)
w_2 = Utils.embedding(B, T, visualize=True, GGH=True, BKZ=False)

# Print the results and the distances between the target point and the found lattice points
print(f"CVP found by R is {w.tolist()}, t-w = {Utils.vector_l2_norm(T-w)}")
print(f"CVP found by B is {w_2.tolist()}, t-w_2 = {Utils.vector_l2_norm(T-w_2)}")

# Note:
# The accuracy of the result using the public basis (B) depends heavily on its Hadamard ratio.
# If the Hadamard ratio of B is very low (< 0.50), the result is likely to be incorrect.
# Even with a good Hadamard ratio, there may still be errors due to the heuristic nature of the method.