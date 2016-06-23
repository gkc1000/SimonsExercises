#!/usr/bin/python

from pyscf import gto

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
print "Number of basis fns=", nao
# overlap, kinetic, nuclear attraction
s = mol.intor('cint1e_ovlp_sph')
t = mol.intor('cint1e_kin_sph')
v = mol.intor('cint1e_nuc_sph')

# The one-electron part of the Hamiltonian is
# the sum of the kinetic and nuclear integrals
h = t + v
    
# 2e integrals (electron repulsion integrals, "eri")
eri = mol.intor('cint2e_sph')

# ERI is stored as [pq, rs] 2D matrix
print "ERI shape=", eri.shape

# Reshape it into a [p,q,r,s] 4D array
eri = eri.reshape([nao,nao,nao,nao])
print "ERI shape=", eri.shape

