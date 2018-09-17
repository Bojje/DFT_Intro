import os
os.chdir("/home/bojje/Dropbox/Fysik og Nanoteknologi/Bachelor/Code/Tutorials/")

from ase import Atoms
from gpaw import GPAW

d = 1.1   # bondlength of hydrogen molecule
a = 5.0   # sidelength of unit cell
c = a / 2
atoms = Atoms('CO',
              positions=[(c - d / 2, c, c),
                         (c + d / 2, c, c)],
              cell=(a, a, a))

calc = GPAW(nbands=5, h=0.2, txt=None)
atoms.set_calculator(calc)

# Start a calculation:
energy = atoms.get_potential_energy()

# Save wave functions:
calc.write('CO.gpw', mode='all')


#from __future__ import print_function
# Setting working directory:
from ase.io import write
from gpaw import restart

basename = 'CO'

# load binary file and get calculator
atoms, calc = restart(basename + '.gpw')

# loop over all wfs and write their cube files
nbands = calc.get_number_of_bands()
for band in range(nbands):
    wf = calc.get_pseudo_wave_function(band=band)
    fname = '{0}_{1}.cube'.format(basename, band)
    print('writing wf', band, 'to file', fname)
    write(fname, atoms, data=wf)
