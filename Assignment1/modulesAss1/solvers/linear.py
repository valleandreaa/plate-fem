#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:09:53 2020

@author: angelo.simone@unipd.it
"""
import sys
import numpy as np
from numpy.linalg import inv
from modulesAss1.plots import *
from modulesAss1.FEM_engine import *
from modulesAss1.funtions import *

# -----------------------------------------------------------------------------

def linear(mesh, BCs, MaterialSets, parameters, alpha):

    print("   . inizializing global arrays")
    K = np.zeros(shape=(mesh.dofs,mesh.dofs)) # allocate memory for global stiffness matrix
    F = np.zeros(shape=(mesh.dofs,1)) # allocate memory for global stiffness matrix

    print("   . generating global stiffness matrix")
    for e in range(mesh.Elements):

        k = stiffness_matrix(e,mesh,MaterialSets,parameters, alpha)
        #print(k)
        # extract system dofs associated with element
        dofs = DofMap(e,mesh)
        # print("Global to local mapping:",dof)
    
        # assemble element matrices into global stiffness matrix
        K = assemble(K,k,dofs)

        # print('dof', dof)





    print("   . applying boundary conditions")
    K, F =BCs.apply(K, F, mesh)

    U = np.matmul(inv(K), F)  # equivalent to inv(K).dot(F)
    print("   . solving KU=F for U")
    return U, K
        