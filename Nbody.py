from environment import *
from amuse.support import io
import sys
import numpy
import numpy as np
from matplotlib import pyplot as plt
from amuse.units import (units, constants)
from amuse.ext.orbital_elements import new_binary_from_orbital_elements
from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.composition_methods import *
from amuse.lab import Particles
from amuse.lab import nbody_system
from amuse.lab import new_powerlaw_mass_distribution
from amuse.plot import scatter
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.couple import bridge
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

BBH0 = io.read_set_from_file(
           'BBH0.txt', 'txt',
           attribute_types = (units.MSun, units.pc, units.kms, units.kms, units.kms, units.pc, units.pc, units.pc),
           attribute_names= ("mass", "radius", "vx", "vy", "vz", "x", "y", "z")
    )

n_BH = 1000 # Number of stellar-mass BHs
Mmin = 1 
Mmax = 100

BBH_nbody = []

i = 0

while i<10:
    print(f"\nTrial {i:d}:")
    try:
        BBH_info = grav(n_BH=n_BH, Mmin=Mmin, Mmax=Mmax)
        
    except IOError:
        print("Encounter some errors and retry:\n")
        continue
    i += 1
    BBH_nbody.append(BBH_info)
    
np.save("Results/Nbody/BBH_nbody", np.array(BBH_nbody))
print("Results saved in Results/Nbody/.")