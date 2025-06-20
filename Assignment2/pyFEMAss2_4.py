"""
@author: andrea.valle.3@studenti.unipd.it
"""
import numpy as np
import sys

# Path modules
sys.path.append("./FEM/Assignment2/modulesAss2")

import meshio
import modulesAss2.BC_engine as BC_engine
import modulesAss2.mesh_engine as mesh_engine
import modulesAss2.FEM_engine as FEM_engine
import modulesAss2.solvers.solver as solver
import modulesAss2.functions as functions
# Dictionary for caracterization of the problem
Procedures = {"solver": {"type": "linear"},
            }


mesh_file = "./meshio/MeshFiles/Piastra_script_100"
mesh = mesh_engine.GMSH(mesh_file)


# Setting the paramenters
MaterialSets = {
    'plane': {'plane deformation': 'plane strain',
            'material behavior': 'isotropic linear elastic',
            'elastic properties': {"Young's modulus": 2.1*pow(10,7), "Poisson's ratio": 0.29},
            'geometric properties': {'thickness': 1},
          # 'stiffness matrix': {'evaluation': 'closed form'}
          'stiffness matrix': {'evaluation': 'numerical integration',
                   'domain': 'line-tri',
                   'rule': 'GaussLegendre',
                   'points': 1}
            },
    }

# Setting the boundary conditions
BCs = BC_engine.BoundaryConditions()     

bc_sym_y = BCs.set("Dirichlet", "sym_y", mesh, [0], [0.0])

# Enforcing displacement at the symmetric nodes
bc_sym_y_enf=[]
for sym in  bc_sym_y:
    Ux, Uy = functions.EnforcingSymDispl(mesh,xCoord=mesh.points[sym[0]][0], yCoord=mesh.points[sym[0]][1])
    sym[3]= Ux
    bc_sym_y_enf.append(sym)

bc_sym_x = BCs.set("Dirichlet", "sym_x", mesh, [1], [ 0.0])
bc_sym_x_enf=[]

for sym in  bc_sym_x:
    Ux, Uy = functions.EnforcingSymDispl(mesh,xCoord=mesh.points[sym[0]][0], yCoord=mesh.points[sym[0]][1])
    sym[3]= Uy
    bc_sym_x_enf.append(sym)

NNodesFrontal= len(BCs.set("Neumann", "frontal", mesh, [0], [1.0]))


bc_frontal = BCs.set("Neumann", "frontal", mesh, [0], [1.0])
bc_frontal = FEM_engine.NodalEquivalentLoads(bc_frontal, mesh)



BCs.data = bc_sym_y_enf + bc_sym_x_enf + bc_frontal


# Solving the problem
U = solver.run(mesh, BCs, MaterialSets, Procedures)



U = U.reshape(len(mesh.points),mesh.d)

stress, VonMises, errorVonMises, VonMisesExact = FEM_engine.triangle_stress(mesh, U, MaterialSets, Procedures)

mesh.point_data = {'Displacements': U.reshape(len(mesh.points),mesh.d)}
err = functions.error(errorVonMises, mesh)
print(err,'error')
mesh.cell_data = { 'VonMises':[VonMises]}

cells = {'triangle': mesh.elements}


meshio.write_points_cells(mesh_file+".vtk", mesh.points, cells, mesh.point_data, mesh.cell_data, binary=False)

