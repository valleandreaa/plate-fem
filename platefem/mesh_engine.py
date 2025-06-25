import numpy as np
import os
import sys
import meshio

def triangle_skewness(p):
    """Return skewness of a triangular element.

    Parameters
    ----------
    p : ndarray (3, 2)
        Coordinates of the triangle vertices.

    Returns
    -------
    float
        Skewness value based on the maximum deviation from the ideal 60 degree
        angle normalised by 60.
    """

    a = np.linalg.norm(p[1] - p[0])
    b = np.linalg.norm(p[2] - p[1])
    c = np.linalg.norm(p[0] - p[2])

    if min(a, b, c) == 0:
        return 1.0

    angles = []
    cosA = (b ** 2 + c ** 2 - a ** 2) / (2 * b * c)
    cosB = (c ** 2 + a ** 2 - b ** 2) / (2 * c * a)
    cosC = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
    for cosv in (cosA, cosB, cosC):
        cosv = np.clip(cosv, -1.0, 1.0)
        angles.append(np.arccos(cosv) * 180.0 / np.pi)

    return max(abs(angles[0] - 60.0), abs(angles[1] - 60.0), abs(angles[2] - 60.0)) / 60.0


def quad_skewness(p):
    """Return skewness of a quadrilateral element.

    Parameters
    ----------
    p : ndarray (4, 2)
        Coordinates of the quad vertices ordered consecutively.

    Returns
    -------
    float
        Skewness value based on the maximum deviation from the ideal 90 degree
        angle normalised by 90.
    """

    angles = []
    for i in range(4):
        v1 = p[i - 1] - p[i]
        v2 = p[(i + 1) % 4] - p[i]
        n1 = np.linalg.norm(v1)
        n2 = np.linalg.norm(v2)
        if n1 == 0 or n2 == 0:
            return 1.0
        cosv = np.dot(v1, v2) / (n1 * n2)
        cosv = np.clip(cosv, -1.0, 1.0)
        angles.append(np.arccos(cosv) * 180.0 / np.pi)

    return max(abs(a - 90.0) for a in angles) / 90.0

def get_key(my_dict, val):
    """
    Function to return key for any value
    """

    for key, value in my_dict.items():
        if val == value[0]:
            return key

    sys.exit()


def _compute_skewness(mesh):
    """Compute skewness for all elements of the mesh."""
    skew = []
    for etype, conn in zip(mesh.elementType, mesh.elements):
        pts = mesh.points[conn]
        if etype == "triangle":
            skew.append(triangle_skewness(pts))
        elif etype == "quad":
            skew.append(quad_skewness(pts))
    return np.array(skew)


def GMSH(mesh_file, element_type="triangle"):

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

    if quad and element_type == "quad":
        meshing = True
        quads = len(mesh.cell_data_dict["gmsh:physical"]["quad"])
        for q in range(quads):
            mesh.elements.append(mesh.cells_dict["quad"][q])
            materialTag = mesh.cell_data_dict["gmsh:physical"]["quad"][q]
            key = get_key(mesh.field_data, materialTag)
            mesh.elementMaterialTag.append(key)
            mesh.elementType.append("quad")
        mesh.NodesElement = 4

    elif quad and element_type != "quad":
        print("quadrilateral elements present but 'element_type' is triangle")
        sys.exit()

    if triangle and element_type == "triangle":
        meshing = True
        triangles = len(mesh.cell_data_dict["gmsh:physical"]["triangle"])
        for t in range(triangles):
            mesh.elements.append(mesh.cells_dict["triangle"][t])
            materialTag = mesh.cell_data_dict["gmsh:physical"]["triangle"][t]

            key = get_key(mesh.field_data, materialTag)
            mesh.elementMaterialTag.append(key)
            mesh.elementType.append("triangle")

    elif triangle and element_type == "quad" and not quad:
        print("Requested quad elements but mesh contains only triangles")
        sys.exit()

    if not meshing:
        print("something went wrong: could not extract mesh data")
        sys.exit()

    mesh.skewness = _compute_skewness(mesh)

    return mesh
