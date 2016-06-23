#!/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#

from pyscf import gto
from pyscf import scf
from pyscf import mcscf
from pyscf import dmrgscf

import os
from pyscf.dmrgscf import settings
'''
Use DMRG program as the solver 
The DMRG program is invoked through system call.  
'''
b = 0.7
mol = gto.M( verbose = 4,
                     atom = [['H', (0, 0, 0)],
                             ['H', (r, 0, 0)],
                             ['H', (2*r, 0, 0)],
                             ['H', (3*r, 0, 0)],
                             ['H', (4*r, 0, 0)],
                             ['H', (5*r, 0, 0)],
                             ['H', (6*r, 0, 0)],
                             ['H', (7*r, 0, 0)],
                             ['H', (8*r, 0, 0)],
                             ['H', (9*r, 0, 0)]],
             basis = 'sto-3g',
)
mf = scf.RHF(mol)
mf.kernel()

mc = dmrgscf.dmrgci.DMRGCI(mf, 8, 8)
mc.mo_coeff = mf.mo_coeff # sets orbitals with which to do DMRG calculation (just HF MO here)
emc = mc.kernel()[0]
print(emc)
