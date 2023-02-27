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

from environment import *

N_BBH = 20
BHdisk, Pinner1, converter1 = make_BHdisk_around_SMBH(SMBH0, 200)

BBH0 = Particles(0)
n = 0

print("Generating initial BBHs:")
for i in range(1, 201):
#     print(n,"/("+str(i)+"/200)")
    ProgressBar(20, n)
    if n >= N_BBH:
        break
    BBH_test = BBH0.copy()
    
    while True:
        flag = 1
        
        BBHi_1mass, BBHi_2mass = new_powerlaw_mass_distribution(2, 1.0|units.MSun, 
                                                          100.0|units.MSun, -2.35)
        np.random.seed()
        a = (np.random.rand(1)[0]*500+500) | units.au
        np.random.seed()
        e = np.random.rand(1)[0]/2
        BBHi = new_binary_from_orbital_elements(BBHi_1mass, BBHi_2mass, 
                                              a, e, G=constants.G)
        #setattr(BBHi, "name", [f"BBH{n+1:d}_1", f"BBH{n+1:d}_2"])

        BBHi.position += BHdisk[BHdisk.name==f'BH{i:d}'].position
        
        # avoid BBHs coming too close
        if len(BBH0)>0: 
            for j in range(0, len(BBH0)):
                if abs(relative(BBHi.center_of_mass(), SMBH0.position).value_in(units.pc)-
                       relative(BBH0[j].position, SMBH0.position).value_in(units.pc))<0.3:
                    flag = 0
                    break 
        if flag==0:
            break
        
        BBHi.velocity += BHdisk[BHdisk.name==f'BH{i:d}'].velocity       
        BBHi.radius = 2*constants.G*BBHi.mass/constants.c**2
        BBH_test.add_particles(BBHi)
        hardness = get_hardness(BBHi, BBH_test)  
        
        # constraint the initial hardness of BBH
        if (hardness>=5) & (hardness<9): 
            BBH0 = BBH_test.copy()
            n += 1
            break
        else:
            BBH_test = BBH0.copy()
            
print(int(len(BBH0)/2), "soft BBHs created.")


from amuse.support import io
io.write_set_to_file(BBH0, 'BBH0.txt', 'txt',
           attribute_types = (units.MSun, units.pc, units.kms, units.kms, units.kms, units.pc, units.pc, units.pc),
           attribute_names= ("mass", "radius", "vx", "vy", "vz", "x", "y", "z")
    )

# Evolve with only SMBH+BBHs
SMBH_ref = SMBH0.copy()
BBH_ref = BBH0.copy()
bodies_ref = Particles(0)
bodies_ref.add_particles(SMBH_ref)
bodies_ref.add_particles(BBH_ref)

BBH_ref = bodies_ref[1:]
SMBH_ref = bodies_ref[0]

gravity = ph4(converter1, number_of_workers=32)
gravity.particles.add_particles(bodies_ref)
channel = gravity.particles.new_channel_to(bodies_ref)

# Stopping condition of collision 
stopping_condition = gravity.stopping_conditions.collision_detection
stopping_condition.enable()
collision_radius_multiplication_factor = 1e5

bodies_ref.radius *= collision_radius_multiplication_factor

# Evolving parameters
end_time = 10 | units.Myr
model_time = 0 | units.Myr
dt = 1e-2 | units.Myr
x = 0

BBH_ref_info = []
print("Evolve with only SMBH+BBHs as reference:")
while(model_time<end_time):
    
    ProgressBar(end_time.value_in(units.Myr), model_time.value_in(units.Myr))
    
    gravity.evolve_model(model_time)
    resolve_collision(stopping_condition, gravity, bodies_ref) # Collisions detection
    channel.copy()
    
    if model_time.value_in(units.Myr)>=x:
        x += 0.5
        
        BBH_ref_info_time = []
        for i in range(N_BBH):
            name1, name2 = BBH_ref[i*2:i*2+2].name
            mass1, mass2 = BBH_ref[i*2:i*2+2].mass.value_in(units.MSun)
            a = BBH_ref[i*2].position - BBH_ref[i*+2+1].position
            hardness = get_hardness(BBH_ref[i*2:i*2+2], bodies_ref)
            a = relative(BBH_ref[i*2].position, BBH_ref[i*2+1].position).value_in(units.au)
            r = relative(BBH_ref[i*2:i*2+2].center_of_mass(), bodies_ref.center_of_mass()).value_in(units.pc)
            v1_v2 = relative(BBH_ref[i*2].velocity, BBH_ref[i*2+1].velocity).value_in(units.km/units.s)
            
            BBH_ref_info_time.append([model_time.value_in(units.Myr), 
                                  name1, name2, mass1, mass2, hardness, a, r, v1_v2])
        
        BBH_ref_info.append(BBH_ref_info_time)
        
    model_time += dt

gravity.stop()
BBH_ref_info = np.array(BBH_ref_info)
np.save("Results/Reference/BBH_ref_info", BBH_ref_info)

print("\nPlot hardness variation...")
plt.figure()
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(BBH_ref_info[:,:,5].astype("float"),axis=1), c="b")
plt.xlabel('Time [Myr]', size=12)
plt.ylabel(r'Average Hardness', size=12)
#plt.legend(bbox_to_anchor=(1.0,1.0))
plt.xlim([0,9.5])
plt.title("Average Hardness of BBHs in Reference")
plt.savefig("Figures/Reference/hardness_ref_average.png", dpi=200, bbox_inches='tight')
plt.show()

for i in range(N_BBH):
    plt.plot(BBH_ref_info[:,i,0].astype('float'), BBH_ref_info[:,i,5].astype('float'))
    plt.xlabel('Time [Myr]')
    plt.ylabel('Hardness')
    plt.savefig(f"Figures/Reference/hardness of BBH_{i:d}.png", dpi=100, bbox_inches='tight')
    plt.close()
