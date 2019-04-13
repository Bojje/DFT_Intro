from contextlib import contextmanager
import os,sys
from ase.io import read, write
from ase.db import connect
from ase.constraints import UnitCellFilter
from ase.optimize.bfgs import BFGS
from gpaw import *
from pathlib import Path
import numpy as np

material = 'V2I4_3653.xyz'
atoms = read('/home/niflheim/s164000/spin_coupl/materials/{0}.xyz'.format(material)) 


#input for relaxation
ecut = 800
# kpoints using the N_k * L approx 20 Å
# kpoints = [5,3,1]
# kpoints using (roughly) the same density as in the c2db
kpoints = [10, 6, 1]

#make sure that periodic boundary conditions and vacuum is correct
atoms.center(vacuum = 8, axis = 2)
atoms.set_pbc([True,True,True])

#setting up a relax command
def relax(atoms, name, calc_mag_state):
    # relax only along xy-axis to opitmize lattice
    sf = UnitCellFilter(atoms,mask=[1,1,0,0,0,0])
    opt = BFGS(sf,trajectory ='{0}_relaxed.traj'.format(name), logfile = '{0}_relaxed.log'.format(name))
    opt.run(fmax=0.05) # until forces < 0.05 eV/atom
    calc.write('relaxed.gpw') #write gpw file for the magnetic relaxations 

calc = GPAW(txt='{0}.txt'.format(name,magstate),
            mode=PW(ecut),
            xc='LDA',
            basis='dzp', #start basis er dzp (nødvendigt for beskrivelse af magnetisme)
            kpts={'size': kpoints, 'gamma': True},
            occupations=FermiDirac(width= 0.05))
# Setting up the parallelspins:
atoms.set_calculator(calc)
# Performing the relaxation
relax(atoms, name, magstate)
