#!/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#
'''
Solve FCI problem with given 1-electron and 2-electron Hamiltonian
'''

import numpy
from pyscf import fci

# This generates a set of random 1e, 2e integrals
numpy.random.seed(12)
norb = 12 # 12 orbitals 
h1 = numpy.random.random((norb,norb))
h2 = numpy.random.random((norb,norb,norb,norb))
# Restore permutation symmetry
h1 = h1 + h1.T
h2 = h2 + h2.transpose(1,0,2,3)
h2 = h2 + h2.transpose(0,1,3,2)
h2 = h2 + h2.transpose(2,3,0,1)

# If you are sure the system ground state is singlet, you can use the spin0 solver.
# The spin0 solver takes this spin symmetry into account to reduce computation cost.
#
cisolver = fci.direct_spin0.FCISolver()
cisolver.verbose = 5
e, fcivec = cisolver.kernel(h1, h2, norb, 8) # 8 particles
print "Energy=", e
