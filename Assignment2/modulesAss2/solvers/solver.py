"""
@author: andrea.valle.3@studenti.unipd.it
"""


def run(mesh, BCs, MaterialSets, Procedures):

    solver_type = Procedures["solver"]["type"]

    from importlib import import_module
    module = import_module(f"{__package__}.{solver_type}")
    solve_fn = getattr(module, solver_type)

    U = solve_fn(mesh, BCs, MaterialSets)

    return U