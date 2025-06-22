"""
@author: andrea.valle.3@studenti.unipd.it
"""


def run(mesh, BCs, MaterialSets, Procedures):

    solver_type = Procedures["solver"]["type"]

    module = __import__(f"platefem.solvers.{solver_type}", fromlist=[solver_type])
    solve_fn = getattr(module, solver_type)

    U = solve_fn(mesh, BCs, MaterialSets)

    return U