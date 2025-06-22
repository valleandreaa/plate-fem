"""Command line interface for building and solving plate problems."""

import meshio

from . import BC_engine, FEM_engine, mesh_engine
from .solvers import solver


def input_nodes():
    """Interactively collect node coordinates from the user."""

    nodes = []
    print("Enter node coordinates 'x y' (type 'done' to finish):")
    while True:
        line = input('node> ').strip()
        if line.lower() == 'done':
            break
        try:
            x, y = map(float, line.split())
        except ValueError:
            print("Invalid input. Provide two numbers or 'done'.")
            continue
        nodes.append((x, y))
    return nodes


def input_lines():
    """Interactively collect line and arc definitions."""

    lines = []
    arcs = []
    print("Define segments with 'n1 n2 name'.")
    print("Define arcs with 'arc n1 center n2 name'.")
    print("Type 'done' when finished.")
    while True:
        line = input('line> ').strip()
        if line.lower() == 'done':
            break
        parts = line.split()
        if not parts:
            continue
        if parts[0].lower() == 'arc':
            if len(parts) != 5:
                print("Format for arc: arc n1 center n2 name")
                continue
            _, n1, center, n2, name = parts
            arcs.append((int(n1), int(center), int(n2), name))
        else:
            if len(parts) != 3:
                print("Format for segment: n1 n2 name")
                continue
            n1, n2, name = parts
            lines.append((int(n1), int(n2), name))
    return lines, arcs


def write_geo(filename, nodes, lines, arcs, mesh_size, surf_name):
    """Write a gmsh .geo file describing the geometry."""
    with open(filename + '.geo', 'w') as f:
        f.write(f"lc = {mesh_size};\n")
        for i, (x, y) in enumerate(nodes, start=1):
            f.write(f"Point({i}) = {{{x}, {y}, 0, lc}};\n")
        line_id = 0
        line_groups = {}
        all_line_ids = []
        for n1, n2, name in lines:
            line_id += 1
            f.write(f"Line({line_id}) = {{{n1}, {n2}}};\n")
            line_groups.setdefault(name, []).append(line_id)
            all_line_ids.append(line_id)
        for n1, c, n2, name in arcs:
            line_id += 1
            f.write(f"Circle({line_id}) = {{{n1}, {c}, {n2}}};\n")
            line_groups.setdefault(name, []).append(line_id)
            all_line_ids.append(line_id)
        if all_line_ids:
            ids = ','.join(str(i) for i in all_line_ids)
            f.write(f"Line Loop(1) = {{{ids}}};\n")
            f.write("Plane Surface(1) = {1};\n")
            f.write(f"Physical Surface('{surf_name}') = {{1}};\n")
        for name, ids in line_groups.items():
            ids_str = ','.join(str(i) for i in ids)
            f.write(f"Physical Curve('{name}') = {{{ids_str}}};\n")


def main():
    """Run the interactive CLI for creating and solving a plate model."""
    nodes = input_nodes()
    if len(nodes) < 3:
        print("Need at least 3 nodes to create a surface")
        return
    lines, arcs = input_lines()
    if not lines and not arcs:
        print("No lines defined")
        return
    mesh_size = float(input("Characteristic length for mesh: "))
    surf_name = input("Name for the surface region: ") or "surface"
    mesh_file = "generated_mesh"
    write_geo(mesh_file, nodes, lines, arcs, mesh_size, surf_name)

    print("Meshing with gmsh...")
    mesh = mesh_engine.GMSH(mesh_file)

    MaterialSets = {
        '1': {
            'plane deformation': input("Plane deformation (plane stress/plane strain): "),
            'material behavior': 'isotropic linear elastic',
            'elastic properties': {"Young's modulus": float(input("Young's modulus: ")),
                                   "Poisson's ratio": float(input("Poisson's ratio: "))},
            'geometric properties': {'thickness': float(input("Thickness: "))},
            'stiffness matrix': {'evaluation': 'closed form'}
        }
    }

    BCs = BC_engine.BoundaryConditions()
    bc_data = []
    boundary_names = {name for _, _, name in lines} | {name for _, _, _, name in arcs}
    for bname in boundary_names:
        btype = input(f"Boundary '{bname}' type (none/Dirichlet/Neumann): ").strip().lower()
        if btype == 'none' or btype == '':
            continue
        if btype == 'dirichlet':
            ux = float(input("  ux value: "))
            uy = float(input("  uy value: "))
            bc_data += BCs.set('Dirichlet', bname, mesh, [0,1], [ux, uy])
        elif btype == 'neumann':
            fy = float(input("  force value in y: "))
            bc_data += BCs.set('Neumann', bname, mesh, [1], [fy])
    BCs.data = bc_data

    Procedures = {"solver": {"type": input("Solver type (linear): ") or "linear"}}

    U = solver.run(mesh, BCs, MaterialSets, Procedures)
    U = U.reshape(len(mesh.points), mesh.dofsNode)

    stress, VonMises, _, _ = FEM_engine.triangle_stress(mesh, U, MaterialSets, Procedures)

    mesh.point_data = {'Displacements': U}
    mesh.cell_data = {'VonMises': [VonMises]}
    cells = {'triangle': mesh.elements}
    vtk_file = mesh_file + '.vtk'
    meshio.write_points_cells(vtk_file, mesh.points, cells, point_data=mesh.point_data, cell_data=mesh.cell_data, binary=False)
    print(f"Results written to {vtk_file}")


if __name__ == '__main__':
    main()