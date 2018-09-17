# Creates: atomization.txt
from __future__ import print_function

from ase import Atoms, Atom
from ase.parallel import paropen as open
from gpaw import GPAW, PW, FermiDirac


a = 10.  # Size of unit cell (Angstrom)
c = a / 2

# Hydrogen atom:
atom = Atoms('H',
             positions=[(c, c, c)],
             magmoms=[0],
             cell=(a, a + 0.0001, a + 0.0002))  # Break cell symmetry

# gpaw calculator:
calc = GPAW(mode=PW(), 
            xc='PBE', 
            hund=True, 
            eigensolver='rmm-diis',  # This solver can parallelize over bands
            occupations=FermiDirac(0.0, fixmagmom=True), 
            txt='H.out',
            )
atom.set_calculator(calc)

e1 = atom.get_potential_energy()
calc.write('H.gpw')

# Hydrogen molecule:
d = 0.74  # Experimental bond length
molecule = Atoms('H2',
                 positions=([c - d / 2, c, c],
                            [c + d / 2, c, c]),
                 cell=(a, a, a))

calc.set(txt='H2.out')
calc.set(hund=False)  # No hund rule for molecules

molecule.set_calculator(calc)
e2 = molecule.get_potential_energy()
calc.write('H2.gpw')

from decimal import Decimal

fd = open('atomization.txt', 'w')
print('  hydrogen atom energy:     %5.2f eV' % e1, file=fd)
print('  hydrogen atom energy:     {0:.4} eV'.format(e1))
print('  hydrogen molecule energy: %5.2f eV' % e2, file=fd)
print('  hydrogen molecule energy: {0:.2} eV'.format(e2))
print('  atomization energy:       %5.2f eV' % (2 * e1 - e2), file=fd)
print('  atomization energy:       {0:.2} eV'.format((2 * e1 - e2)))
fd.close()

#Running the optimization script:
# Creates: optimization.txt
#from __future__ import print_function

from gpaw import restart
from ase.parallel import paropen as open
from ase.optimize import QuasiNewton
from ase.constraints import FixAtoms

molecule, calc = restart('H2.gpw', txt='H2-relaxed.txt')
molecule.set_constraint(FixAtoms(mask=[0, 1]))

e2 = molecule.get_potential_energy()
d0 = molecule.get_distance(0, 1)

fd = open('optimization.txt', 'w')
print('experimental bond length:', file=fd)
print('hydrogen molecule energy: %5.2f eV' % e2, file=fd)
print('bondlength              : %5.2f Ang' % d0, file=fd)

# Find the theoretical bond length:
relax = QuasiNewton(molecule, logfile='qn.log')
relax.run(fmax=0.05)

e2 = molecule.get_potential_energy()
d0 = molecule.get_distance(0, 1)

print(file=fd)
print('PBE energy minimum:', file=fd)
print('hydrogen molecule energy: %5.2f eV' % e2, file=fd)
print('bondlength              : %5.2f Ang' % d0, file=fd)
fd.close()
