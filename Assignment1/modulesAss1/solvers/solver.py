"""
@author: andrea.valle.3@studenti.unipd.it
"""


def run(mesh, BCs, MaterialSets, parameters, alpha):

    solverType = parameters["solver"]["type"]
    
    exec("from modulesAss1.solvers."+solverType+" import "+solverType)

    U, K = eval (solverType+"(mesh, BCs, MaterialSets, parameters, alpha)")

    
    return U, K

