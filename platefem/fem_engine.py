"""
@author: andrea.valle.3@studenti.unipd.it
"""
import numpy as np
import sys

import math
from numpy.linalg import multi_dot

from .FEM_utils import *
from .FEM_utils import quadrature_rule, shape_functions
from .functions import *
from .material_data import *
def stiffness_matrix(e, mesh, MaterialSets):
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
    elType = mesh.elementType[e]
    elMat = mesh.elementMaterialTag[e]
    key = MaterialSets[str(elMat)]

    elementPoints = mesh.points[mesh.elements[e]][:,:mesh.d]

    nodesInElement = len(elementPoints) # or len(mesh.elements[e])
    elementDofs = nodesInElement * mesh.dofsNode
    matrix = np.zeros((elementDofs, elementDofs))

    if mesh.d == 1:
        strainComponents = 1
    elif mesh.d == 2:
        strainComponents = 3
    elif mesh.d == 3:
        strainComponents = 6
    else:
        print("Error, spatial dimension different from 1, 2, 3: mesh.d=", mesh.d)
        sys.exit()

    D = np.zeros((strainComponents, strainComponents))

    if key['stiffness matrix']["evaluation"] == 'numerical integration':
        evaluation = "numerical integration"
        domain = key['stiffness matrix']["domain"]
        rule = key['stiffness matrix']["rule"]
        points = key['stiffness matrix']["points"]
    elif key['stiffness matrix']["evaluation"] == 'closed form':
        evaluation = "closed form"
        domain = None
        rule = None
        points = None
    else:
        raise TypeError("Keyword not recognized for stiffness matrix evaluation:\n",
                        key['stiffness matrix']["evaluation"])

    if evaluation == 'closed form':


        if elType == "triangle":
            E = key['elastic properties']["Young's modulus"]
            nu = key['elastic properties']["Poisson's ratio"]
            thickness = key['geometric properties']["thickness"]

            # plane deformation: plane stress or plane strain?
            if key["plane deformation"] == 'plane stress':
                D = E / (1 - nu ** 2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2.]])
            elif key["plane deformation"] == 'plane strain':
                D = E / (1 + nu) / (1 - 2 * nu) * np.array(
                    [[1 - nu, nu, 0], [nu, 1 - nu, 0], [0, 0, (1 - 2 * nu) / 2.]])
            else:
                print("stiffness_matrix(): closed form: Cannot define D!")
                sys.exit()

            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
            n3 = mesh.elements[e][2]
            x1 = mesh.points[n1, 0]
            y1 = mesh.points[n1, 1]
            x2 = mesh.points[n2, 0]
            y2 = mesh.points[n2, 1]
            x3 = mesh.points[n3, 0]
            y3 = mesh.points[n3, 1]

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


        if elType == "triangle":
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

        for ip in range(len(points)):
            # Evaluate determinant of Jacobian matrix and derivatives of shape functions

            detJ, dN, N, B = shape_functions(mesh.elements[e], points[ip], mesh.NodesElement, mesh.d, mesh.points)

            matrix += multi_dot([B.T, D, B]) * detJ/2 * weights[0]*thickness


        return matrix

    else:

        print("you should not be here----")
        sys.exit()


def stress_fem(e, mesh, MaterialSets, U):
    """
    Retrieving the stresses
    e: number of element
    mesh: innstance of the mesh
    MaterialSets: conditions of the model
    U: displacements
    """

    elMat = mesh.elementMaterialTag[e]
    key = MaterialSets[str(elMat)]
    elType = mesh.elementType[e]
    evaluation = key['stiffness matrix']["evaluation"]
    dofs = DofMap(e, mesh)

    if evaluation == 'closed form':

        if elType == "triangle":
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

            n1 = mesh.elements[e][0]
            n2 = mesh.elements[e][1]
            n3 = mesh.elements[e][2]
            x1 = mesh.points[n1, 0]
            y1 = mesh.points[n1, 1]
            x2 = mesh.points[n2, 0]
            y2 = mesh.points[n2, 1]
            x3 = mesh.points[n3, 0]
            y3 = mesh.points[n3, 1]

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
            u = np.zeros(6)

            u[0]= U[mesh.elements[e][0]][0]
            u[1]= U[mesh.elements[e][0]][1]
            u[2]= U[mesh.elements[e][1]][0]
            u[3]= U[mesh.elements[e][1]][1]
            u[4]= U[mesh.elements[e][2]][0]
            u[5]= U[mesh.elements[e][2]][1]

            Strain =multi_dot([B, u])

            Stress= multi_dot([D, Strain])


            return Stress

    if evaluation == 'numerical integration':
        if elType == "triangle":
            E = key['elastic properties']["Young's modulus"]
            nu = key['elastic properties']["Poisson's ratio"]
            thickness = key['geometric properties']["thickness"]

            # plane deformation: plane stress or plane strain?
            if key["plane deformation"] == 'plane stress':
                D = E / (1 - nu ** 2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2.]])
            elif key["plane deformation"] == 'plane strain':
                D = E / (1 + nu) / (1 - 2 * nu) * np.array(
                    [[1 - nu, nu, 0], [nu, 1 - nu, 0], [0, 0, (1 - 2 * nu) / 2.]])
            else:
                print("stiffness_matrix(): closed form: Cannot define D!")
                sys.exit()

            d=mesh.d
            NodesElement=mesh.NodesElement
            pointsEl = mesh.points
            N = np.zeros((d, NodesElement))
            dN = np.zeros((d, NodesElement))
            J = np.zeros((d, d))
            elementPoints = np.zeros((3, 2))
            element_coordinates=mesh.elements[e]
            n1 = element_coordinates[0]
            n2 = element_coordinates[1]
            n3 = element_coordinates[2]
            x1 = pointsEl[n1, 0]
            y1 = pointsEl[n1, 1]
            x2 = pointsEl[n2, 0]
            y2 = pointsEl[n2, 1]
            x3 = pointsEl[n3, 0]
            y3 = pointsEl[n3, 1]
            elementPoints[0][0] = x1
            elementPoints[0][1] = y1
            elementPoints[1][0] = x2
            elementPoints[1][1] = y2
            elementPoints[2][0] = x3
            elementPoints[2][1] = y3

            csi = 1/3
            eta = 1/3
            N1 = csi
            N2 = eta
            N3 = 1 - csi - eta

            N = np.array([[N1, 0, N2, 0, N3, 0], [0, N1, 0, N2, 0, N3]])
            dN = np.array([[1, 0, -1],
                           [0, 1, -1]])

            J = np.dot(dN, elementPoints)
            detJ = determinant(J)

            invJ = inverse(J)

            dN = np.dot(invJ, dN).reshape((d, NodesElement))
            B = np.array([[dN[0][0], 0, dN[0][1], 0, dN[0][2], 0],
                          [0, dN[1][0], 0, dN[1][1], 0, dN[1][2]],
                          [dN[1][0], dN[0][0], dN[1][1], dN[0][1], dN[1][2], dN[0][2]]])
            u = np.zeros(6)

            u[0] = U[mesh.elements[e][0]][0]
            u[1] = U[mesh.elements[e][0]][1]
            u[2] = U[mesh.elements[e][1]][0]
            u[3] = U[mesh.elements[e][1]][1]
            u[4] = U[mesh.elements[e][2]][0]
            u[5] = U[mesh.elements[e][2]][1]

            Strain = multi_dot([B, u])

            Stress = multi_dot([D, Strain])

        return Stress


    else:

        print("you should not be here-----")
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
    nodesInElement = len(mesh.elements[e])
    nodesElement = mesh.NodesElement
    dofs = [mesh.elements[e][i] * dofsNode + j for i in range(nodesInElement) for j in range(dofsNode)]
    return dofs

def triangle_stress(mesh, U, MaterialSets, Procedures):
    """
    Triangle
    """
    points = mesh.points
    Stress = []
    VonMises =[]
    VonMisesExactlist =[]
    error = 0
    for e in range(len(mesh.elements)):
        stress = stress_fem(e,mesh,MaterialSets,U)
        Stress.append(stress)

        Ux, Uy, Sigma_xx, Sigma_yy, Sigma_xy, VonMisesExact = ExactCalc(mesh,e)
        VonMises.append((math.sqrt(Stress[e][0] ** 2 + Stress[e][1] ** 2 -
                         Stress[e][0] * Stress[e][1] +
                         3 * Stress[e][2]**2)))
        VonMisesExactlist.append(VonMisesExact)
    VonMises = np.array(VonMises)
    VonMisesExactlist =np.array(VonMisesExactlist)
    errorVonMises=VonMisesExactlist -VonMises
    # print('VonMisesExactlist',VonMisesExactlist)
    # print('VonMises',VonMises)





    return Stress, VonMises, errorVonMises, VonMisesExactlist

def NodalEquivalentLoads(edge, mesh):
    load=100
    thickness=1
    EdgeDist=abs(mesh.points[edge[0][0]][1]-mesh.points[edge[1][0]][1])/(len(edge)-1)

    F_nodal=0.5*load*EdgeDist*thickness

    for i in range(len(edge)):
        if mesh.points[edge[i][0]][1]==0 or mesh.points[edge[i][0]][1]==10:

            edge[i][3] = F_nodal
        else:
            edge[i][3] = 2*F_nodal
    return edge