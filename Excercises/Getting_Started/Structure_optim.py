import os
os.chdir('/home/bojje/Dropbox/Fysik og Nanoteknologi/Bachelor/Code/Excercises/Getting_Started/')

from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
#
# Eksemple udregning
#
system = Atoms('H2', positions=[[0.0, 0.0, 0.0],
                                [0.0, 0.0, 1.0]])
calc = EMT()

system.set_calculator(calc)

opt = QuasiNewton(system, trajectory='h2.emt.traj')

opt.run(fmax=0.05)

# Lad os prøve for H_2O

system = Atoms('H2O', positions=[[0.0, 0.0, 0.0],
                                 [0.0, 0.0, 1.0],
                                 [0.0, 1.0, 0.0]])
calc = EMT()

system.set_calculator(calc)

opt = QuasiNewton(system, trajectory='h2O.emt.traj')

opt.run(fmax=0.05)

#
# Brug af gpaw i stedet for
#
from ase import Atoms
from gpaw import GPAW
from ase.optimize import QuasiNewton

system = Atoms('H2O', positions=[[0.0, 0.0, 0.0],
                                 [0.0, 0.0, 1.0],
                                 [0.0, 1.0, 0.0]])
system.set_cell((6.0, 6.0, 6.0))
system.center()

calc = GPAW()
system.set_calculator(calc)

opt = QuasiNewton(system, trajectory='h2oGPAW.emt.traj')

opt.run(fmax=0.05)
