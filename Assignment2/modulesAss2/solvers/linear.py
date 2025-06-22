#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:09:53 2020

@author: angelo.simone@unipd.it
"""

import numpy as np
from numpy.linalg import inv

from ..FEM_engine import *

# -----------------------------------------------------------------------------

def linear(mesh, BCs, MaterialSets):

    print("   . inizializing global arrays")
    systemDofs = mesh.dofsNode*len(mesh.points) # total number of dofs
    K = np.zeros(shape=(systemDofs,systemDofs)) # allocate memory for global stiffness matrix
    F = np.zeros(shape=(systemDofs,1)) # allocate memory for global stiffness matrix

    print("   . generating global stiffness matrix")

    for e in range(len(mesh.elements)):

        k = stiffness_matrix(e,mesh,MaterialSets)
        #print(k)

        # extract system dofs associated with element
        dofs = DofMap(e,mesh)
        # print("Global to local mapping:",dof)
    
        # assemble element matrices into global stiffness matrix
        K = assemble(K,k,dofs)

        # print('dof', dof)




    print("   . applying boundary conditions")
    K, F =BCs.apply(K, F, mesh)

    print("   . solving KU=F for U")
    U = np.matmul(inv(K),F) # equivalent to inv(K).dot(F)

    return U
        