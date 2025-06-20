import numpy as np
import sys

def geometric_properties(MaterialSets, materialNumber, elType):
    """
    Retrieving geometric properties
    MaterialSets: Dictionary with caracteristics of the problem to solve
    materialNumber: specify the integer of the material
    elType: speficy the type of element evalueted
    """

    key = MaterialSets[str(materialNumber)]


    if key['element'] == 'bar' and elType == 'bar':
        area = key['geometric properties']['area']
        return area

    else:
        print("geometric_properties(): keyword not allowed")
        sys.exit()


def elastic_properties(MaterialSets, materialNumber, elType):
    """
    Retrieving geometric propertieselastic properties
    MaterialSets: Dictionary with caracteristics of the problem to solve
    materialNumber: specify the integer of the material
    elType: speficy the type of element evalueted
    """


    key = MaterialSets[str(materialNumber)]

    if key['element'] == 'bar' and elType == 'bar':
        Young = key['elastic properties']["Young's modulus"]
        return Young

    else:
        print("elastic_properties(): keyword not allowed")
        sys.exit()

def stiffness_matrix_evaluation_info(MaterialSets, materialNumber, elType):
    """
    A function that returns information about the stiffness matrix evaluation.
    key: entiring the materialsets
    evaluation: type of finite elements
    domain: type of integration
    rule: rule of integration
    points: number of points integration
    """
    key = MaterialSets[str(materialNumber)]['stiffness matrix']

    if key["evaluation"] == 'numerical integration':
        evaluation = "numerical integration"
        domain = key["domain"]
        rule = key["rule"]
        points = key["points"]

    elif key["evaluation"] == 'closed form':
        evaluation = "closed form"
        domain = None
        rule = None
        points = None

    else:
        print("stiffness_matrix_properties(): I don't know how this stiffness matrix is created. What does", key,
              " mean?")
        sys.exit()

    return evaluation, domain, rule, points