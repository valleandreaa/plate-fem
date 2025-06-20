"""
@author: andrea.valle.3@studenti.unipd.it
"""


def run(mesh, BCs, MaterialSets, Procedures):

    solverType = Procedures["solver"]["type"]
    
    exec("from modulesAss2.solvers."+solverType+" import "+solverType)

    U = eval (solverType+"(mesh, BCs, MaterialSets)")
    

    
    return U

