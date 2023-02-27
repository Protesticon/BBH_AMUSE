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
            BBH0, 'BBH0.txt', 'txt',
           attribute_types = (units.MSun, None, units.pc, units.kms, units.kms, units.kms, units.pc, units.pc, units.pc),
           attribute_names= ("mass", "name", "radius", "vx", "vy", "vz", "x", "y", "z")
    )