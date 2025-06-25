"""
@author: andrea.valle.3@studenti.unipd.it
"""
import numpy as np
import sys


def stiffness_matrix_evaluation_info(MaterialSets, materialNumber, elType):
    """
    A function that returns information about the stiffness matrix evaluation.
    key: entiring the materialsets
    evaluation: type of finite elements
    domain: type of integration
    rule: rule of integration
    points: number of points integration
    """

    key = MaterialSets[str(int(materialNumber))]['stiffness matrix']

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