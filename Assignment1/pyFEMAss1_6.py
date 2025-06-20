"""
@author: andrea.valle.3@studenti.unipd.it
"""


import numpy as np
from numpy.linalg import inv
from modulesAss1.FEM_engine import *
from modulesAss1.mesh_engine import *
from modulesAss1.material_data import *
from modulesAss1.BC_engine import *
from modulesAss1.solvers.solver import *
from modulesAss1.plots import *
from modulesAss1.funtions import *
# Setting the scenarios
# Number of Elements
ListElements=[4, 8, 16, 32, 64]
# Values of the parameters alpha
alphaList=[0.5,1,2,4]
# Setting how to solve the problem, one of the following
NumericalIntegration=True
for alpha in alphaList:
    # Lists for keeping errors and size of discratization
    err = []
    hSize = []
    for NElements in ListElements:
        # Calling instances
        mesh = Mesh()
        BCs = BoundaryConditions()
        # Dimensions of the problem to solve
        mesh.d = 1
        mesh.NodesElement = 2
        # Number of elements
        mesh.Elements = NElements
        # Dictionary for caracterization of the problem
        parameters = {"strain components": mesh.d,
                      "solver": {"type": "linear"},
                      }

        mesh.elements = np.zeros((NElements, 4))
        # Definition of the total number of nodes
        mesh.Nodes = mesh.Elements + 1
        mesh.NodesElement = 2
        mesh.dofsNode = 1
        # Total number of degree of freedom
        mesh.dofs = mesh.dofsNode * mesh.Nodes

        # Array with the coordinates of the nodes, not equispatial this time
        for j in range(NElements):
            mesh.elements[j]=(int(2), int(1), int(j), int(j+1))
        mesh.elements = mesh.elements.astype(int)
        nodes = NElements+1
        factor = 0
        for x in range(nodes): factor += x
        h = 10 / factor
        mesh.coordinates = np.zeros((nodes,))
        point = 0
        for x in range(nodes):
            point += x * h
            mesh.coordinates[x] = point


        if NumericalIntegration == True:
            MaterialSets = {
                '1': {'element': 'bar',
                      'elastic properties': {"Young's modulus": 1},
                      'geometric properties': {'area': 1},
                      'stiffness matrix': {'evaluation': 'numerical integration',
                                           'domain': 'line',
                                           'rule': 'GaussLegendre',
                                           'points': 2}
                      },
            }
        # Defining the boundary conditions as integer
        Neumann=0
        Dirichlet=1
        # Getting the node corrisponding to x=10
        LastNode= mesh.Nodes -1
        # List for the Boundary conditions
        BCs.data=[
                 [LastNode, Dirichlet, 0, -0.000045399],
                 [0, Neumann  , 0, -1.0],
                 ]
        # Solving the problem
        U, K= run(mesh, BCs, MaterialSets, parameters, alpha)

        # Plotting the 2-nodes bar solution
        if mesh.NodesElement ==2 and alpha ==1:
            PlotDisplacement(x=mesh.coordinates, y=U, i=NElements, Xlabel='x', Ylabel='U(x)', NameFile='./Img/Assignment_1_3_2Nodes.png')

        # Evalueting the error from the discretization
        epsilon_h = 0.5 * multi_dot([U.T, K, U])

        err, hSize  = ErrorFunction(err, hSize,h, alpha=alpha, P=1, epsilon_h=epsilon_h)

    PlotErrorAlpha(x=[5, 9, 17, 33, 65], y=err, label=str(alpha), Xlabel='Number of nodes', Ylabel='Error',
              NameFile='./Img/errorRefinement.png')


    if alpha ==1:PlotError(x=hSize, y=err, Xlabel='h Size', Ylabel='Error',
                           NameFile='./Img/error_alpha1Refinement.png')



