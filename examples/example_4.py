
import numpy as np
import sys
import meshio
from platefem import bc_engine, mesh_engine, fem_engine, functions
from platefem.solvers import solver

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
BCs = bc_engine.BoundaryConditions() 

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
bc_frontal = fem_engine.NodalEquivalentLoads(bc_frontal, mesh)



BCs.data = bc_sym_y_enf + bc_sym_x_enf + bc_frontal


# Solving the problem
U = solver.run(mesh, BCs, MaterialSets, Procedures)



U = U.reshape(len(mesh.points),mesh.d)

stress, VonMises, errorVonMises, VonMisesExact = fem_engine.triangle_stress(mesh, U, MaterialSets, Procedures)

mesh.point_data = {'Displacements': U.reshape(len(mesh.points),mesh.d)}
err = functions.error(errorVonMises, mesh)
print(err,'error')
mesh.cell_data = { 'VonMises':[VonMises]}

cells = {'triangle': mesh.elements}


meshio.write_points_cells(mesh_file+".vtk", mesh.points, cells, mesh.point_data, mesh.cell_data, binary=False)

