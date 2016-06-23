import numpy
from pyscf import gto, scf, fci, mp, cc

for r in numpy.arange(0.5, 3.5, .1):
    mol = gto.M(
        atom = [['H', (0, 0, 0)],
                ['H', (r, 0, 0)],
                ['H', (2*r, 0, 0)],
                ['H', (3*r, 0, 0)],
                ['H', (4*r, 0, 0)],
                ['H', (5*r, 0, 0)]],
        basis = '6-31g',
        verbose = 0
    )
    myhf = scf.RHF(mol)
    myhf.kernel()

    # FCI calc
    myci = fci.FCI(mol, myhf.mo_coeff)
    e, civec = myci.kernel()

    # MP2 calc
    mymp = mp.MP2(myhf)
    emp2, t2 = mymp.kernel() # T2 is the set of MP2 amplitudes (amplitudes in the excitation operator)

    # CCSD calc
    mycc = cc.CCSD(myhf)
    ecc, t1, t2 = mycc.kernel() # T1, T2 are the CC amplitudes (amplitudes in the excitation operator)

    print('r = %f, EFCI = %g, ECC = %g, EMP = %g' % (r, e+mol.energy_nuc(), mycc.ecc+myhf.energy_tot(), mymp.emp2 + myhf.energy_tot()))
