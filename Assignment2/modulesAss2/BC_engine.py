import sys
import numpy as np
from numpy import unique

# -----------------------------------------------------------------------------

class BoundaryConditions:
    """
    Contains information about the boundary conditions
    data: list containing the boundary condiutions
    mesh: class of the mesh engine
    """

    def __init__(self):
        self.data = None

    # -----------------------------------------------------------------------------

    def apply(self, K, F, mesh):

        '''
        Definition of the boundary conditions
        K: Global Stiffness matrix
        F: vector of forces
        dof: global degree of freedom
        '''


        for bc in self.data:  # this select a row at a time in data

            if bc[1] == 'Neumann':  # Neumann, 0
                dof = bc[0] * mesh.dofsNode + bc[2]
                F[dof] += bc[3]

            elif bc[1] == 'Dirichlet':  # Dirichlet, 1
                dof = bc[0] * mesh.dofsNode + bc[2]
                F -= K[:, [dof]] * bc[3]
                K[:,dof] = np.zeros(len(K))
                K[dof,:] = np.zeros(len(K))
                K[dof, dof] = 1  # set diagonal to 1
                F[dof] = bc[3]  # enforce value
            else:
                print("BC not recognized")
                sys.exit(0)

        return K, F

    def set(self, BCtype, name, mesh, dof, value):
        """
        Retrieving data from the mesh instance
        BCtype: type of boundary conditions
        name:
        mesh: mesh instance
        dof: number degree of freedom
        value: value of the boundary condition
        """
        tag = mesh.field_data[name][0]
        dim = mesh.field_data[name][1]

        if dim == 0:
            # array containing indices of elements on the boundary
            on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["vertex"] == tag)[0]
            # array containing indices of nodes on the boundary
            nodes = unique(mesh.cells_dict["vertex"][on_boundary])
        elif dim == 1:
            on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["line"] == tag)[0]
            nodes = unique(mesh.cells_dict["line"][on_boundary])
        else:
            print("Error in BC_engine.set(): dimension not coded")

        listBCs = []

        for n in nodes:
            for d in range(len(dof)):
                listBCs.append([n, BCtype, dof[d], value[d]])

        return listBCs
