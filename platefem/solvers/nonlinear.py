import numpy as np
from numpy.linalg import inv

from ..FEM_engine import stiffness_matrix, assemble, DofMap


def nonlinear(mesh, BCs, MaterialSets, max_iter=20, tol=1e-6):
    """Simple Newton-Raphson nonlinear solver."""
    system_dofs = mesh.dofsNode * len(mesh.points)
    U = np.zeros((system_dofs, 1))
    for iteration in range(max_iter):
        K = np.zeros((system_dofs, system_dofs))
        F = np.zeros((system_dofs, 1))
        for e in range(len(mesh.elements)):
            k = stiffness_matrix(e, mesh, MaterialSets)
            dofs = DofMap(e, mesh)
            K = assemble(K, k, dofs)
        K, F = BCs.apply(K, F, mesh)
        R = F - K @ U
        if np.linalg.norm(R) < tol:
            print(f"   . converged in {iteration} iterations")
            break
        dU = inv(K) @ R
        U += dU
    else:
        print("   . warning: nonlinear solver did not converge")
    return U