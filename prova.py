from flint import *
from GGH_crypto import Utils

mat = Utils.load_matrix_from_file("totry.txt")

Utils.BKZ_reduction(mat, block=60, precision=90, pruned=True)