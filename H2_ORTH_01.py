from pyscf import gto, scf, ao2mo, molden, localizer

r = 0.7 # Angstrom

mol = gto.M( 
    atom = [['H', (0, 0, 0)],
            ['H', (r, 0, 0)],
            ['H', (2*r, 0, 0)],
            ['H', (3*r, 0, 0)]],
    basis = 'sto-3g',
    verbose = 0
)
nao = mol.nao_nr()
print "Number of basis fns =", nao

mf = scf.RHF(mol)
mf.kernel() 
print('E(HF) = %.12f' % (mf.energy_tot()))

# Use Molden to visualize these files
filename_mo = 'H4.molorbs.molden'
filename_boys = 'H4.boysorbs.molden'

with open( filename_mo, 'w' ) as thefile:
    molden.header( mol, thefile )
    molden.orbital_coeff( mol, thefile, mf.mo_coeff )
    print("Molecular orbitals saved in", filename_mo)

    # try these different options
    # tolocalize = [0,1,2,3] # localize all orbitals together
    tolocalize = [0,1] # localize only occupied
    # tolocalize = [2,3] # localize only virtual
    loc  = localizer.localizer( mol, mf.mo_coeff[:,tolocalize], 'boys' )

    new_coeff = loc.optimize()
    loc.dump_molden( filename_boys, new_coeff )
    print("Boys localized orbitals saved in", filename_boys)

    print "Columns are the localized orbital coefficents"
    print new_coeff
