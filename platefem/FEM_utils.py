"""
@author: andrea.valle.3@studenti.unipd.it
"""
import numpy as np
import sys

# -----------------------------------------------------------------------------

def quadrature_rule(domain, rule, points):
    """
    Define quadrature rules.
    The quadrature rule is returned as a tuple (of tuples).
    rule: type of quadrature
    domain: type of element
    points: number of points and then position of the points
    weight: weight of the points
    """

    # TODO: the nested-if approach below is not pretty but it works. Better solution?
    # TODO: more robust check on error?
    # TODO: increase number of digits

    if rule == "GaussLegendre":
        if domain == "line":
            if points == 1:
                weights = (2.0)
                points = (0.0)
            elif points == 2:
                weights = (1.0, 1.0)
                points = (-0.5773502691896257,
                          +0.5773502691896257)
            elif points == 3:
                weights = (0.5555555555555556,
                           0.8888888888888888,
                           0.5555555555555556)
                points = (-0.7745966692414834,
                          0.0,
                          +0.7745966692414834)
        if domain == "line-tri":
            if points == 1:
                weights = [0.5]
                points = (1 / 3,
                          1 / 3)


    return weights, points



# -----------------------------------------------------------------------------

def shape_functions(element_coordinates, points, NodesElement, d,pointsEl):
    '''
    Evaluate derivatives of shape functions at the integration point location.
    element_coordinates:
    N: array of shape functions
    dN: derivative of the shape functions
    J: jacobian matrix
    point (xi): coordinate of the element in local reference

    '''


    # Allocation memory

    d=int(d)
    N = np.zeros((d,NodesElement))
    dN = np.zeros((d,NodesElement))
    J = np.zeros((d,d))
    elementPoints = np.zeros((3,2))



    if NodesElement==3 and d ==2:
        n1 = element_coordinates[0]
        n2 = element_coordinates[1]
        n3 = element_coordinates[2]
        x1 = pointsEl[n1, 0]
        y1 = pointsEl[n1, 1]
        x2 = pointsEl[n2, 0]
        y2 = pointsEl[n2, 1]
        x3 = pointsEl[n3, 0]
        y3 = pointsEl[n3, 1]
        elementPoints[0][0]=x1
        elementPoints[0][1]=y1
        elementPoints[1][0]=x2
        elementPoints[1][1]=y2
        elementPoints[2][0]=x3
        elementPoints[2][1]=y3
        csi = points
        eta = points
        N1 = csi
        N2 = eta
        N3 = 1 - csi - eta

        N = np.array([[N1, 0, N2, 0, N3, 0], [0, N1, 0, N2, 0, N3]])
        dN = np.array([[1, 0, -1],
                       [0, 1, -1]])


        J = np.dot(dN,elementPoints)

        detJ = determinant(J)

        invJ = inverse(J)

        dN = np.dot(invJ, dN).reshape((d, NodesElement))
        B=np.array([[dN[0][0], 0,       dN[0][1], 0,        dN[0][2],0],
                    [ 0, dN[1][0],       0, dN[1][1],        0,dN[1][2]],
                    [dN[1][0], dN[0][0], dN[1][1], dN[0][1], dN[1][2], dN[0][2]]])


    if detJ < 0:
        print("Error: negative Jacobian")
        sys.exit()



    return detJ, dN, N, B


def determinant(matrix):
    """
    Determinat of a matrix
    matrix: selected matrix
    """
    det = np.linalg.det(matrix)
    return det

def inverse(matrix):
    """
    Inverse of a matrix
    matrix: selected matrix
    """
    inv = np.linalg.inv(matrix)
    return inv