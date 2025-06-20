import numpy as np
from numpy.linalg import multi_dot
from numpy.linalg import inv

# Young's modulus
E=1
# Area
A=1
# Integration points. It is available 1, 2 or 3 integration method
points = 3
if points == 1:
    w = np.array([2])  # weight
    x = np.array([0])  # position points
if points == 2:
    w = np.array([1, 1])  # weight
    x = np.multiply(3 ** (-0.5), [-1, 1])  # position points

if points == 3:
    w = np.array([5 / 9, 8 / 9, 5 / 9])  # weight
    x = np.multiply((3 / 5) ** (0.5), [-1, 0, 1])  # position points
# inizialization shape function vector
N = np.zeros((1, 3))
# inizialization derivative shape function vector
dN = np.zeros((1, 3))
# inizialization jacobian matrix
J = np.zeros((3, 3))
# inizialization stiffness shape function matrix
K = np.zeros((3, 3))
# inizialization coordiante of the nodes vector
el_coord = np.zeros((1, 3))
for e in range(points):
    # shape functions 3-nodes element
    N[0][0] = 0.5 * x[e] * (x[e] - 1)
    N[0][1] = 1 - x[e] ** 2
    N[0][2] = 0.5 * x[e] * (x[e] + 1)

    # derivation of the shape function
    dN[0][0] = x[e] - 0.5
    dN[0][1] = -2 * x[e]
    dN[0][2] = x[e] + 0.5

    # definition nodes' position
    el_coord[0][0] = -1
    el_coord[0][1] = 0
    el_coord[0][2] = +1
    # jacobian matrix
    J = np.dot(dN, el_coord.reshape(3, 1))
    # determinant of the jacobian
    detJ = np.linalg.det(J)
    # inverse of the jacobian
    invJ = np.linalg.inv(J)
    # redifinition of the jacobian matrix
    dN = np.dot(invJ, dN).reshape((1, 3))
    # derivation of the shape function vector
    B = np.array([dN[0][0], dN[0][1], dN[0][2]]).reshape(1, 3)
    # assembly global stiffness matrix
    K += multi_dot([B.T, B]) * E * A * detJ * w[e]
print('Stiffness matrix 3-nodes element:', K)