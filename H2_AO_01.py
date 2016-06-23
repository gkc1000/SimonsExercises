#!/usr/bin/python

from pyscf import gto, scf, fci

r = 0.7 # Angstrom

# This defines a molecule object, which
# holds the geometry, basis, nuclear repulsion energy, etc.
# See also the examples in PATH_TO_PYSCF/examples/gto/
mol = gto.M( 
    atom = [['H', (0, 0, 0)],
            ['H', (r, 0, 0)]],
    basis = 'sto-3g',
    verbose = 0
)
# number of basis functions
nao = mol.nao_nr()
print "Number of basis fns =", nao

# set up a Hartree-Fock calculation object
myhf = scf.RHF(mol)
myhf.kernel() # executes the calculation
              # returns various information (see scf/hf.py)
              # as well as stores it in the RHF object
print('E(HF) = %.12f' % (myhf.energy_tot()))
#
# Function fci.FCI creates an FCI solver based on the given orbitals and the
# num. electrons and spin of the given mol object
#
cisolver = fci.FCI(mol, myhf.mo_coeff)
print('E(FCI) = %.12f' % (cisolver.kernel()[0] + mol.energy_nuc()))
              
