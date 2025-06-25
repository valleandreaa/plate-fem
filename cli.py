"""Run a plate analysis using parameters from a YAML file."""

import sys
import yaml
import meshio
import matplotlib.pyplot as plt

from Assignment2.modulesAss2 import BC_engine, FEM_engine, mesh_engine
from Assignment2.modulesAss2.solvers import solver
from cli import write_geo


def plot_setup(mesh, BCs):
    """Visualize mesh nodes and applied boundary conditions."""
    x = mesh.points[:, 0]
    y = mesh.points[:, 1]
    tris = mesh.elements
    plt.triplot(x, y, tris, color="lightgray")

    dir_nodes = {bc[0] for bc in BCs.data if bc[1] == 'Dirichlet'}
    neu_nodes = {bc[0] for bc in BCs.data if bc[1] == 'Neumann'}

    if dir_nodes:
        plt.scatter(x[list(dir_nodes)], y[list(dir_nodes)], color='red', label='Dirichlet')
    if neu_nodes:
        plt.scatter(x[list(neu_nodes)], y[list(neu_nodes)], color='green', label='Neumann')
        for bc in BCs.data:
            if bc[1] == 'Neumann':
                n = bc[0]
                dof = bc[2]
                scale = 0.05 * (y.max() - y.min())
                if dof == 1:
                    plt.arrow(x[n], y[n], 0, bc[3]*scale, color='blue', head_width=scale*0.2)
                elif dof == 0:
                    plt.arrow(x[n], y[n], bc[3]*scale, 0, color='blue', head_width=scale*0.2)
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Plate configuration')
    plt.axis('equal')
    plt.show()


def run_from_yaml(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)

    nodes = cfg['nodes']
    lines = cfg.get('lines', [])
    arcs = cfg.get('arcs', [])
    mesh_conf = cfg['mesh']
    mesh_file = mesh_conf.get('file', 'generated_mesh')
    write_geo(mesh_file, nodes, lines, arcs, mesh_conf['size'], mesh_conf.get('surface_name', 'surface'), mesh_conf.get('element_type', 'triangle'))

    mesh = mesh_engine.GMSH(mesh_file, mesh_conf.get('element_type', 'triangle'))

    MaterialSets = cfg['materials']

    BCs = BC_engine.BoundaryConditions()
    bc_data = []
    for name, binfo in cfg.get('boundaries', {}).items():
        btype = binfo['type'].lower()
        if btype == 'dirichlet':
            bc_data += BCs.set('Dirichlet', name, mesh, [0, 1], binfo['values'])
        elif btype == 'neumann':
            bc_data += BCs.set('Neumann', name, mesh, [1], [binfo['value']])
    BCs.data = bc_data

    Procedures = cfg.get('procedures', {'solver': {'type': 'linear'}})

    U = solver.run(mesh, BCs, MaterialSets, Procedures)
    U = U.reshape(len(mesh.points), mesh.dofsNode)

    stress, VonMises, _, _ = FEM_engine.triangle_stress(mesh, U, MaterialSets, Procedures)

    mesh.point_data = {'Displacements': U}
    mesh.cell_data = {'VonMises': [VonMises]}
    cells = {'triangle': mesh.elements}
    vtk_file = mesh_file + '.vtk'
    meshio.write_points_cells(vtk_file, mesh.points, cells, point_data=mesh.point_data, cell_data=mesh.cell_data, binary=False)
    print(f"Results written to {vtk_file}")

    plot_setup(mesh, BCs)


def main():
    if len(sys.argv) != 2:
        print("Usage: python yaml_cli.py <config.yaml>")
        return
    run_from_yaml(sys.argv[1])


if __name__ == '__main__':
    main()