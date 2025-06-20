import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import multi_dot
from numpy.linalg import inv


def StiffnessMatrixTri(w, csi,eta,ni):
    E = 1  # Young's modulus
    D = (E / (1 - ni ** 2)) * np.array([[1, ni, 0],  # matrix of constant elastic
                                        [ni, 1, 0],
                                        [0, 0, (1 - ni) / 2]])



    N = np.zeros((2, 6))  # initialization shape function matrix
    dN = np.zeros((2, 3))  # initialization derivation shape function matrix
    J = np.zeros((2, 3))  # initialization jacobian matrix
    K = np.zeros((3, 3))  # initialization stiffness matrix
    el_coord = np.zeros((2, 3))  # initialization coordinate vector

    N1 = csi
    N2 = eta
    N3 = 1 - csi - eta
    N = np.array([[N1, 0, N2, 0, N3, 0], [0, N1, 0, N2, 0, N3]])
    dN = np.array([[1, 0, -1], [0, 1, -1]])
    el_coord = np.array([[0, 1, 0], [0, 0, 1]])

    J = np.dot(dN, el_coord.T)
    detJ = np.linalg.det(J)
    invJ = np.linalg.inv(J)
    dN = np.dot(invJ, dN).reshape((2, 3))
    # Contract from DN
    B = np.array([[dN[0][0], 0, dN[0][1], 0, dN[0][2], 0],
                  [0, dN[1][0], 0, dN[1][1], 0, dN[1][2]],
                  [dN[1][0], dN[0][0], dN[1][1], dN[0][1], dN[1][2], dN[0][2]]])

    K = multi_dot([B.T, D, B]) * detJ * w
    return K, w, csi, eta, ni
# Show that the result is independent of the integration point location

K, w, csi, eta, ni= StiffnessMatrixTri(w=1/2, csi=1/3, eta=1/3, ni=0)
print('weight:',w,'\n',
      'Integration point coordinates:\n',
      'Csi:', csi,'\n',
      'Eta:', eta,'\n',
      'Ni:', ni,'\n',
      'Siffness Matrix:', K,'\n')

K, w, csi, eta, ni= StiffnessMatrixTri(w=1/2, csi=1/4, eta=1/4, ni=0)
print('weight:',w,'\n',
      'Integration point coordinates:\n',
      'Csi:', csi,'\n',
      'Eta:', eta,'\n',
      'Ni:', ni, '\n',
      'Siffness Matrix:', K, '\n')

# Is this result valid for any value of the Poisson’s ratio?
K, w, csi, eta, ni= StiffnessMatrixTri(w=1/2, csi=1/3, eta=1/3, ni=0.3)
print('weight:',w,'\n',
      'Integration point coordinates:\n',
      'Csi:', csi,'\n',
      'Eta:', eta,'\n',
      'Ni:', ni,'\n',
      'Siffness Matrix:', K,'\n')

K, w, csi, eta, ni= StiffnessMatrixTri(w=1/2, csi=1/4, eta=1/4, ni=0.3)
print('weight:',w,'\n',
      'Integration point coordinates:\n',
      'Csi:', csi,'\n',
      'Eta:', eta,'\n',
      'Ni:', ni, '\n',
      'Siffness Matrix:', K, '\n')

# Form the right-hand-side vector due to self-weight using one-point Gauss integration.
# Is it also independent of the integration point location? Yes
def ForceTri(K):
    U=np.array([1,0,0,0,0,0])
    F = multi_dot([K,U])


ForceTri(K)
