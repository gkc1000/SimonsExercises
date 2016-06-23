#!/usr/bin/python

from pyscf import gto, scf, ao2mo

r = 0.7 # Angstrom

# This defines a molecule object, which
# holds the geometry, basis, nuclear repulsion energy, etc.
# See also the examples in PATH_TO_PYSCF/examples/gto/
mol = gto.M( 
    atom = [['H', (0, 0, 0)],
            ['H', (r, 0, 0)]],
    basis = 'cc-pVDZ',
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

# There are many different ways to carry out the transformation
# from one set of orbitals to another. PATH_TO_PYSCF/examples/ao2mo shows
# some different ones.

# This transforms the AO integrals (associated with the mol object)
# to the MO coefficients obtained after the HF calc
eri_mo = ao2mo.kernel(mol, myhf.mo_coeff)

# Here is another more manual way.
# 
# c = myhf.mo_coeff # These define the MOs to transform to
# h1e = reduce(numpy.dot, (c.T, myhf.get_hcore(), c)) # this multiplies
#                                                     # the 1e part of the Hamiltonian
#                                                     # in the AO basis by the MO coeffients (c)
# eri = ao2mo.incore.full(myhf._eri, c) # this transforms from the eris in the AO basis (myhf._eri)
#                                       # to eris in the new basis defined by c

# Note the MO integrals are stored in a matrix
# (p>=q, r>=s)
# that relies on 4-fold permutation symmetry
# of the integrals, since
# (pq|rs) = (qp|rs) = (pq|sr) = (qp|sr)
print "ERI shape=", eri_mo.shape

# Create the full integral array (without perm. symmetry)
eri_mo = ao2mo.restore(1, eri_mo, nao)
print "ERI shape=", eri_mo.shape
