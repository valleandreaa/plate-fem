import numpy as np
import os
import sys
import meshio

def get_key(my_dict, val):
    """
    Function to return key for any value
    """

    for key, value in my_dict.items():
        if val == value[0]:
            return key

    sys.exit()


def GMSH(mesh_file):

    os.system("gmsh -2 " + mesh_file + ".geo" + " -o " + mesh_file + ".msh")

    # create a mesh object
    mesh = meshio.read(mesh_file + ".msh")

    pyFEM_MeshAttributes = ["d", "dofsNode", "elements", "elementMaterialTag", "elementType", "points"]

    for attribute in pyFEM_MeshAttributes:
        if attribute in dir(mesh):
            if attribute == "points":
                pass
            else:
                print("Error: meshio already contains the attribute", attribute)
                print("       ...do something!")
                sys.exit()


    mesh.d = 2
    mesh.dofsNode = 2
    mesh.NodesElement = 3
    mesh.elements = []
    mesh.elementMaterialTag = []
    mesh.elementType = []
    meshing = False

    quad = False
    try:
        dummy = mesh.cell_data_dict['gmsh:physical']['quad']
        quad = True
    except KeyError:
        # print("No quadrilateral elements in mesh")
        pass

    triangle = False
    try:
        dummy = mesh.cell_data_dict['gmsh:physical']['triangle']
        triangle = True
    except KeyError:
        # print("No triangular elements in mesh")
        pass

    if quad:
        print("quadrilateral elements not implemented")
        sys.exit()

    if triangle:
        meshing = True
        triangles = len(mesh.cell_data_dict["gmsh:physical"]["triangle"])
        for t in range(triangles):
            mesh.elements.append(mesh.cells_dict["triangle"][t])
            materialTag = mesh.cell_data_dict["gmsh:physical"]["triangle"][t]

            key = get_key(mesh.field_data, materialTag)
            mesh.elementMaterialTag.append(key)
            mesh.elementType.append("triangle")

    if not meshing:
        print("something went wrong: could not extract mesh data")
        sys.exit()

    return mesh
