
import numpy as np
from numpy.linalg import multi_dot
from modulesAss1.FEM_utils import quadrature_rule
from modulesAss1.FEM_utils import shape_functions
from modulesAss1.material_data import *

def stiffness_matrix(e, mesh, MaterialSets, parameters,alpha):
    """
    Generate the stiffness matrix for a spring element
    elType: Definition of the element type
    elMat: material element
    nodesInElement: Evaluation of number of nudes in each element
    stiffness_matrix_evaluation_info: info about resolution methode
    NodesElement: number of nodes for each element
    dofsNode: degree of freedom for each element
    elementsDofs: degree of freefom for each element
    matrix: stifness matrix for each element
    """
    elType = mesh.elementType(e)

    elMat = int(mesh.elements[e][1])


    nodesInElement = mesh.nodesInElement(e)

    evaluation, domain, rule, points = stiffness_matrix_evaluation_info(MaterialSets, elMat, elType)

    strainSize = parameters['strain components']

    D = np.zeros((strainSize, strainSize))
    K = np.zeros((strainSize, strainSize))
    elementDofs = mesh.NodesElement * mesh.dofsNode
    matrix = np.zeros((elementDofs, elementDofs))

    if evaluation == 'closed form':

        if elType == "bar":
            Young = elastic_properties(MaterialSets,elMat, elType)
            area =geometric_properties(MaterialSets,elMat, elType)
            n1 = nodesInElement[0]
            n2 = nodesInElement[1]
            length = mesh.coordinates[n2] - mesh.coordinates[n1]
            if length < 0.0:
                print("stiffness_matrix(): closed form: Oh dear, bar", e, "has a negative length!")
                sys.exit()
            return (Young * area *length *alpha**2* np.array([(1/3,  1/6), (1/6, 1/3)]))+ (Young * area / length)* np.array([(1, - 1), (- 1, 1)])


        elif elType == "triangle":
            key = MaterialSets[str(elMat)]
            E = key['elastic properties']["Young's modulus"]
            nu = key['elastic properties']["Poisson's ratio"]
            thickness = key['geometric properties']["thickness"]

            if key["plane deformation"] == 'plane stress':
                D = E / (1 - nu ** 2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2.]])
            elif key["plane deformation"] == 'plane strain':
                D = E / (1 + nu) / (1 - 2 * nu) * np.array(
                    [[1 - nu, nu, 0], [nu, 1 - nu, 0], [0, 0, (1 - 2 * nu) / 2.]])
            else:
                print("stiffness_matrix(): closed form: Cannot define D!")
                sys.exit()

            n1 = nodesInElement[0]
            n2 = nodesInElement[1]
            n3 = nodesInElement[2]
            x1 = mesh.coordinates[n1, 0]
            y1 = mesh.coordinates[n1, 1]
            x2 = mesh.coordinates[n2, 0]
            y2 = mesh.coordinates[n2, 1]
            x3 = mesh.coordinates[n3, 0]
            y3 = mesh.coordinates[n3, 1]

            y23 = y2 - y3
            y31 = y3 - y1
            y12 = y1 - y2
            x32 = x3 - x2
            x13 = x1 - x3
            x21 = x2 - x1
            x32 = x3 - x2
            y13 = y1 - y3
            x23 = x2 - x3

            area = 0.5 * ((x2 * y3 - x3 * y2) + (x3 * y1 - x1 * y3) + (x1 * y2 - x2 * y1))
            detJ = x13 * y23 - y13 * x23

            B = 1 / detJ * np.array([[y23, 0, y31, 0, y12, 0],
                                     [0, x32, 0, x13, 0, x21],
                                     [x32, y23, x13, y31, x21, y12]])

            matrix = thickness * area * multi_dot([B.T, D, B])

            return matrix
        else:
            print("stiffness_matrix(): closed form: I don't know elType =", elType)
            sys.exit()

    elif evaluation == 'numerical integration':
        weights, points = quadrature_rule(domain, rule, points)

        if elType == "bar":

            Young = elastic_properties(MaterialSets,elMat, elType)
            D[0][0] = Young
            E= 1
            A = 1
            K[0][0] = E*A*alpha**2
            area =geometric_properties(MaterialSets,elMat, elType)
            VolumeFactor = area


        else:
            print("stiffness_matrix(): numerical integration: I don't know elType =", elType)
            sys.exit()

        for ip in range(len(points)):
            detJ,dN, N,B= shape_functions(mesh.element_coordinates(e), points[ip], mesh.NodesElement, mesh.d)

            matrix += multi_dot([B.T, D, B]) * VolumeFactor* detJ*weights[ip]+multi_dot([N.T, N])*K * VolumeFactor * detJ* weights[ip]

        return matrix

    else:

        print("you should not be here")
        sys.exit()

def assemble(K, k, dofs):
    """
    Assemble a stiffness matrix into the global stiffness matrix
    K: global stiffness matrix
    k: local stiffness matrix
    dofs: degree of freedome
    """

    K[np.ix_(dofs, dofs)] += k
    return K


def DofMap(e, mesh):
    """
    Generation of the global stiffness matrix
    """

    dofsNode = mesh.dofsNode
    nodesInElement = mesh.nodesInElement(e)
    nodesElement = mesh.NodesElement

    dofs = [nodesInElement[i] * dofsNode + j for i in range(nodesElement) for j in range(dofsNode)]
    return dofs