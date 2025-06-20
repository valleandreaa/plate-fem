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
ClosedForm =False
NumericalIntegration=True
for alpha in alphaList:
    # Lists for keeping errors and size of discratization
    err=[]
    hSize=[]
    for NElements in ListElements:
        # Calling instances
        mesh=Mesh()
        BCs = BoundaryConditions()
        # Dimensions of the problem to solve
        mesh.d = 1
        # Number of elements
        mesh.Elements = NElements

        # Number nodes for each element (2 or 3)
        mesh.NodesElement = 2 #or3
        # Total of degree of freedom per node
        mesh.dofsNode = 1
        # Dictionary for caracterization of the problem
        parameters = {"strain components": mesh.d,
                      "solver": {"type": "linear"},
                      }

        # Definition of the total number of nodes considering either the 2 or 3 nodes bar
        if mesh.NodesElement ==2:mesh.Nodes = mesh.Elements + 1
        if mesh.NodesElement ==3:mesh.Nodes =int(mesh.Elements*2 + 1)

        # Total number of degree of freedom
        mesh.dofs = mesh.dofsNode * mesh.Nodes
        # Considering a number of nodes equal to 2, here is built the array of the elements
        if mesh.NodesElement == 2:
            # Initialisation array for data elements
            mesh.elements = np.zeros((NElements, 4))
            for j in range(NElements):
                mesh.elements[j]=(int(2), int(1), int(j), int(j+1))
            mesh.elements = mesh.elements.astype(int)
            # Size of the elements
            h =  10/NElements
            # Array with the coordinates of the nodes
            mesh.coordinates = np.linspace(0.0, 10.0, NElements+1)

        # Considering a number of nodes equal to 3, here is built the array of the elements
        if mesh.NodesElement == 3:
            # Initialisation array for data elements
            mesh.elements = np.zeros((NElements, 5))
            l=0
            for j in range(int(NElements)):
                mesh.elements[j] = (int(2), int(1), int(l), int(l + 1), int(l + 2))
                l+=2
            mesh.elements = mesh.elements.astype(int)
            # Size of the elements
            h = 10 / int(NElements*1.5+ 1)
            # Array with the coordinates of the nodes
            mesh.coordinates = np.linspace(0.0, 10.0,int(NElements*1.5+ 1) )

        print("Setting the boundary conditions")

        if ClosedForm==True:
            MaterialSets = {
            '1': {'element': 'bar',
                  'elastic properties': {"Young's modulus": 1},
                  'geometric properties': {'area': 1},
                  'stiffness matrix': {'evaluation': 'closed form'}
                  },
            }

        if NumericalIntegration==True:
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
        # Plotting the 3-nodes bar solution
        if mesh.elements[0][1]==1 and mesh.NodesElement ==3 and alpha ==1:
            PlotDisplacement(x=mesh.coordinates, y=U, i=NElements, Xlabel='x', Ylabel='U(x)', NameFile='./Img/Assignment_1_3_3Nodes.png')
        # Evalueting the error from the discretization
        epsilon_h = 0.5 * multi_dot([U.T, K, U])

        # Getting the error and the size of the elements
        err, hSize  = ErrorFunction(err, hSize,h, alpha=alpha, P=1, epsilon_h=epsilon_h)

    # Plotting the comparison the errors considering differents alpha
    PlotErrorAlpha(x=hSize, y=err, label=str(alpha), Xlabel='h Size', Ylabel='Error',
              NameFile='./Img/error.png')
    # Plotting the errors just for alpha equal to 1
    if alpha ==1:PlotError(x=hSize, y=err, Xlabel='h Size', Ylabel='Error',
                           NameFile='./Img/error_alpha1.png')



