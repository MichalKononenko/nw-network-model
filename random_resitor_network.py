# Random resitor network - 2 point calculation
# Jeremy Smith
# Northwestern University
# Version 1.0

import numpy as np
import turtle as tu


#   N2----N3
#   | \    |
#   |  \   |
#   |   \  |
#   |    \ |
#   N1----N4

# Test case network calculator

r1 = 102.0  # resistor between 1,2 2,3 3,4 1,4
r2 = 5.1    # resistor between 2,4

# Conductance matrix (adjacency)
c_ij = np.array([[0, 1/r1, 0, 1/r1], [1/r1, 0, 1/r1, 1/r2], [0, 1/r1, 0, 1/r1], [1/r1, 1/r2, 1/r1, 0]])

# Conductance matrix (degree)
c_i = np.sum(c_ij, axis=1)

# Laplacian matrix
Lmatrix = c_i*np.identity(len(c_i)) - c_ij

# Find eigenvalues and eigenvectors
val, vec = np.linalg.eigh(Lmatrix)

# Resistances between nodes
R12 = 0
R13 = 0
R24 = 0

for i in range(len(val)):
	if abs(val[i]) < 1e-10:
		continue
	R12 += (1/val[i])*(vec[0][i] - vec[1][i])**2
	R13 += (1/val[i])*(vec[0][i] - vec[2][i])**2
	R24 += (1/val[i])*(vec[1][i] - vec[3][i])**2

print R12
print R13
print R24







