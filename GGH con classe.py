from GGH_crypto import Utils, GGHCryptosystem
from flint import fmpz_mat
from decimal import Decimal
import time
import numpy as np

dimension = 200

GGH_object = GGHCryptosystem(dimension = dimension, debug=True)
GGH_object.encrypt()


