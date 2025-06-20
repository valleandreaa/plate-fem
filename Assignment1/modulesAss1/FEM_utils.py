import numpy as np
import sys

def quadrature_rule(domain, rule, points):
    """
    Define quadrature rules.
    The quadrature rule is returned as a tuple (of tuples).
    rule: type of quadrature
    domain: type of element
    points: number of points and then position of the points
    weight: weight of the points
    """

    if rule == "GaussLegendre":
        if domain == "line":
            if points == 1:
                weights = [2.0]
                points = (0.0)
                points=[0.0]
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



            # else:
            #     print("quadrature_rule(): I don't know any", rule, "rule on", domain, "with", points, "points.")
            #     sys.exit()

        else:
            print("quadrature_rule(): I don't know any", rule, "rule on", domain)
            sys.exit()
    else:
        print("quadrature_rule(): I don't know any", rule, "rule")
        sys.exit()

    return weights, points



# -----------------------------------------------------------------------------

def shape_functions(element_coordinates, integrationPoint, NodesElement,d):
    '''
    Evaluate derivatives of shape functions at the integration point location.
    element_coordinates:
    N: array of shape functions
    dN: derivative of the shape functions
    J: jacobian matrix
    point (xi): coordinate of the element in local reference
    '''


    # Allocation memory
    N = np.zeros((d,NodesElement))
    dN = np.zeros((d,NodesElement))
    J = np.zeros((d,d))
    J_N = np.zeros((d,d))
    elementPoints = element_coordinates.reshape(NodesElement, d)

    if NodesElement==2:
        N[0][0] = -0.5 * integrationPoint + 0.5
        N[0][1] = +0.5 * integrationPoint + 0.5

        dN[0][0] = -0.5
        dN[0][1] = +0.5

        J_N = np.dot(N, elementPoints)
        J = np.dot(dN,elementPoints)

        detJ = determinant(J)
        detJ_N = determinant(J_N)
        invJ = inverse(J)
        dN = np.dot(invJ, dN).reshape((d, NodesElement))
        B = np.array([dN[0][0], dN[0][1]]).reshape(1, 2)

    if NodesElement==3:
        N[0][0] = -0.5 * integrationPoint + 0.5*integrationPoint**2
        N[0][1] = 1 - integrationPoint**2
        N[0][2] = +0.5 * integrationPoint + 0.5*integrationPoint**2


        dN[0][0] = -0.5 +integrationPoint
        dN[0][1] = -2*integrationPoint
        dN[0][2] = +0.5+integrationPoint
        J_N = np.dot(N, elementPoints)
        J = np.dot(dN,elementPoints)

        detJ = determinant(J)
        invJ = inverse(J)

        dN = np.dot(invJ, dN).reshape((d, NodesElement))
        B = np.array([dN[0][0], dN[0][1], dN[0][2] ]).reshape(1, 3)

    if detJ < 0:
        print("Error: negative Jacobian")
        sys.exit()

    return detJ,  dN, N, B


def determinant(matrix):
    det = np.linalg.det(matrix)
    return det

def inverse(matrix):
    inv = np.linalg.inv(matrix)
    return inv