import numpy as np
import sys

class BoundaryConditions:
    """
    Contains information about the boundary conditions
    data: list containing the boundary condiutions
    mesh: class of the mesh engine
    """

    def __init__(self):
        self.data = None


    def apply(self, K, F, mesh):
        '''
        Definition of the boundary conditions
        K: Global Stiffness matrix
        F: vector of forces
        dof: global degree of freedom
        '''
        for bc in self.data:
            if bc[1] == 0:  # Neumann
                dof = bc[0] * mesh.dofsNode + bc[2]
                F[dof] += bc[3]
            elif bc[1] == 1:  # Dirichlet
                dof = bc[0] * mesh.dofsNode + bc[2]
                F -= K[:, [dof]] * bc[3]
                K[:, dof] = np.zeros(mesh.dofs)
                K[dof, :] = np.zeros(mesh.dofs)
                K[dof, dof] = 1
                F[dof] = bc[3]
            else:
                print("BC not recognized")
                sys.exit()

        return K, F