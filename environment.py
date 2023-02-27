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

# Progress bar
def ProgressBar(total, progress):
    barLength, status = 20, ""
    progress = float(progress) / float(total)
    if progress >= 1.:
        progress, status = 1, "\r\n"
    block = int(round(barLength * progress))
    text = "\r[{}] {:.0f}% {}".format(
        "#" * block + "-" * (barLength - block), round(progress * 100, 0),
        status)
    sys.stdout.write(text)
    sys.stdout.flush()
    
def orbital_period(Mtot, a):
    return (((4 * numpy.pi**2) * a**3)/(constants.G * Mtot)).sqrt()

# Define collision detection
def merge_two_stars(bodies, particles_in_encounter):
    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()
    d = (particles_in_encounter[0].position - particles_in_encounter[1].position)
    v = (particles_in_encounter[0].velocity - particles_in_encounter[1].velocity)
    print("Two BHs", particles_in_encounter.name, particles_in_encounter.mass.in_(units.MSun),
          "collided with d=", d.length().in_(units.au))
    
    particles_in_encounter[np.argmax(particles_in_encounter.mass)].mass = particles_in_encounter.total_mass()
    particles_in_encounter[np.argmax(particles_in_encounter.mass)].position = com_pos
    particles_in_encounter[np.argmax(particles_in_encounter.mass)].velocity = com_vel
    particles_in_encounter[np.argmax(particles_in_encounter.mass)].radius = particles_in_encounter.radius.sum()
    bodies.remove_particle(particles_in_encounter[np.argmin(particles_in_encounter.mass)])
    
def resolve_collision(collision_detection, gravity, bodies):
    if collision_detection.is_set():
        print("Well, we have an actual collision between two or more BHs.")
        print("This happened at time=", gravity.model_time.in_(units.Myr))
        for ci in range(len(collision_detection.particles(0))): 
            encountering_particles = Particles(particles=[collision_detection.particles(0)[ci],
                                                          collision_detection.particles(1)[ci]])
            colliding_BHs = encountering_particles.get_intersecting_subset_in(bodies)
            merge_two_stars(bodies, colliding_BHs)
            bodies.synchronize_to(gravity.particles)

            
            
# Binaries detection
def get_binaries(particles,hardness=0.5,G=constants.G):
    n=len(particles)
    
    # define a fixed average kinetic energy
    average_Ek = 1.70 | units.km**2/units.s**2
    
    max_mass=particles.mass.amax()
    limitE=hardness*average_Ek
    
    # sort to speed up loop
    a=numpy.argsort(particles.x.number)
    a=np.delete(a,particles[a].name=="SMBH")
    n=len(a)

    binaries=[]

    for i in range(0, n-1):
        j=i+1
        while j<n and j<i+10 and (particles.x[a[j]]-particles.x[a[i]])<2*G*max_mass/limitE:
            r2=(particles.x[a[j]]-particles.x[a[i]])**2+ \
               (particles.y[a[j]]-particles.y[a[i]])**2+ \
               (particles.z[a[j]]-particles.z[a[i]])**2 
            v2=(particles.vx[a[j]]-particles.vx[a[i]])**2+ \
               (particles.vy[a[j]]-particles.vy[a[i]])**2+ \
               (particles.vz[a[j]]-particles.vz[a[i]])**2 
            r=r2**0.5
            eb=G*(particles.mass[a[i]]+particles.mass[a[j]])/r-0.5*v2
            if eb > limitE:
                binary=particles[[a[i],a[j]]].copy()
                binary.hardness=eb/average_Ek
                com_pos = binary.center_of_mass()
                binary_momentum = binary.velocity * binary.mass.reshape((-1,1))
                L_b = (binary.position-com_pos).cross(binary_momentum).sum(axis=0)
                L_b = L_b.in_(units.pc*units.MSun*units.km/units.s)
                binary.thetaz = L_b[2] / L_b.length() # angle between the angular momentum and z-axis
                binaries.append(binary)
            j+=1  

    return binaries, average_Ek

# Get hardness of a binary in the system
def get_hardness(BBH, particles, G=constants.G):
    
    n=len(particles)
    total_Ek=(0.5*particles.mass*(particles.vx**2+particles.vy**2+particles.vz**2)).sum()
    average_Ek=total_Ek/particles.mass.sum()
    
    r2=(BBH.x[0]-BBH.x[1])**2+ \
       (BBH.y[0]-BBH.y[1])**2+ \
       (BBH.z[0]-BBH.z[1])**2 
    v2=(BBH.vx[0]-BBH.vx[1])**2+ \
       (BBH.vy[0]-BBH.vy[1])**2+ \
       (BBH.vz[0]-BBH.vz[1])**2 
    r=r2**0.5
    eb=G*(BBH.mass[0]+BBH.mass[1])/r-0.5*v2
    average_Ek = 1.70 | (units.km**2/units.s**2)
    hardness = eb/average_Ek
    
    return hardness

# Get relative distance or velocity 
def relative(property1, property2):
    return ((property1 - property2)**2).sum()**0.5

# Start with the SMBH in Sgr A*
SMBH = Particles(1) 
SMBH.name = "SMBH"
SMBH.mass = 4.154e+6 | units.MSun
SMBH.position = (0, 0, 0) | units.pc
SMBH.velocity = (0, 0, 0) | units.kms
SMBH.radius = 2*constants.G*SMBH.mass/constants.c**2
SMBH0 = SMBH.copy()
print("SMBH mass =", SMBH0.mass)
r = 10 # pc (the outermost radius)


# Gas density (rho~rho0*(r/r0)^y, where y(r<r0)=1.2 and y(r>r0)=1.75)
r0 = 0.22 | units.pc
def rho_0(Mdisk, r0=r0):
    return Mdisk/(4*np.pi*r0**3/1.8-4*np.pi*r0**3/1.25+4*np.pi*r0**1.75*(r|units.pc)**1.25/1.25)

# Make BH disk
def make_BHdisk_around_SMBH(SMBH=SMBH0, n_BH=1000, Mmin=1.0, Mmax=100.0):
    R = r | units.pc
    Ndisk = n_BH
    Rin = 0.1
    Rout = 1
    Pinner1 = orbital_period(SMBH.mass, Rin*R)
    converter1 = nbody_system.nbody_to_si(SMBH.mass.sum(), R)
    masses = new_powerlaw_mass_distribution(Ndisk, Mmin|units.MSun, Mmax|units.MSun, -2.35) # BH masses in powerlaw
    BHdisk = ProtoPlanetaryDisk(Ndisk,
                              convert_nbody=converter1,
                              Rmin=Rin,
                              Rmax=Rout,
                              q_out=10.0,
                              discfraction=masses.sum()/SMBH.mass).result
    BHdisk.mass = masses
    BHdisk.name = np.arange(1,n_BH+1).astype("str")
    BHdisk.radius = 2*constants.G*BHdisk.mass/constants.c**2
    BHdisk.remove_attribute_from_store('u')

    for i in range(n_BH):
        BHdisk[i].name = "BH"+BHdisk[i].name
    return BHdisk, Pinner1, converter1

# Make gas disk
def make_gasdisk_around_SMBH(SMBH=SMBH0, Mdisk=1, Ndisk=1000):
    R = r | units.pc
    Ndisk = Ndisk
    Rin = 0.01
    Rout = 1
    Pinner2 = orbital_period(SMBH.mass, Rin*R)
    Mdisk = Mdisk * SMBH.mass
    converter2 = nbody_system.nbody_to_si(SMBH.mass.sum(), R)
    
    gasdisk = ProtoPlanetaryDisk(Ndisk,
                              convert_nbody=converter2,
                              Rmin=Rin,
                              Rmax=Rout,
                              q_out=10.0,
                              discfraction=Mdisk/SMBH.mass).result
    gasdisk.name = "gas disk"

    # gas cloud mass varies with orbital radius (Mtot=Mdisk)
    for i in range(Ndisk):
        r1 = relative(gasdisk[i].position, SMBH.position)
        if r1 < r0:
            gasdisk.mass[i] = (r1/r0)**(-1.2) | units.MSun
        else:
            gasdisk.mass[i] = (r1/r0)**(-1.75) | units.MSun

    gasdisk.mass *= Mdisk/gasdisk.mass.sum()
    
    # gas density assigned as rho0
    rho = rho_0(Mdisk, r0)
    gasdisk.radius = (gasdisk.mass/(4*rho))**(1./3.)
    return gasdisk, Pinner2, converter2


# Evolve with pure Nbody
def grav(n_BH=1000, Mmin=1.0, Mmax=100.0):
    print("Running Nbody (n_BH="+str(n_BH)+", Mmin="+str(round(Mmin,2))+"Msun, Mmax="+str(round(Mmax,2))+"Msun)......")

    BHdisk, Pinner1, converter1 = make_BHdisk_around_SMBH(SMBH0, n_BH, Mmin, Mmax)

    bodies = Particles(0)
    bodies.add_particles(SMBH0)
    bodies.add_particles(BHdisk)
    bodies.add_particles(BBH0)

    SMBH = bodies[bodies.name=="SMBH"][0]
    BHdisk = bodies[1:n_BH+1]
    BBH  = bodies[n_BH+1:]
    
    gravity = ph4(converter1, number_of_workers=32)
    gravity.particles.add_particles(bodies)
    channel = gravity.particles.new_channel_to(bodies)

    # Stopping condition of collision 
    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()
    collision_radius_multiplication_factor = 1e5

    bodies.radius *= collision_radius_multiplication_factor

    # Evolving parameters
    end_time = 10 | units.Myr
    model_time = 0 | units.Myr
    dt = 1e-2 | units.Myr
    x = 0

    BBH_info = []

    # Evolve
    while(model_time<end_time):

        ProgressBar(end_time.value_in(units.Myr), model_time.value_in(units.Myr))

        gravity.evolve_model(model_time)
        
        # collision detection
        resolve_collision(stopping_condition, gravity, bodies)
        channel.copy()

        if model_time.value_in(units.Myr)>=x:
            x += 0.5

            BBH_info_time = []
            for i in range(N_BBH):
                name1, name2 = BBH[i*2:i*2+2].name
                mass1, mass2 = BBH[i*2:i*2+2].mass.value_in(units.MSun)
                a = BBH[i*2].position - BBH[i*+2+1].position
                hardness = get_hardness(BBH[i*2:i*2+2], bodies)
                a = relative(BBH[i*2].position, BBH[i*2+1].position).value_in(units.au)
                r = relative(BBH[i*2:i*2+2].center_of_mass(), SMBH.position).value_in(units.pc)
                v1_v2 = relative(BBH[i*2].velocity, BBH[i*2+1].velocity).value_in(units.km/units.s)

                BBH_info_time.append([model_time.value_in(units.Myr), 
                                      name1, name2, mass1, mass2, hardness, a, r, v1_v2])

            BBH_info.append(BBH_info_time)

        model_time += dt
        
    gravity.stop()
    return np.array(BBH_info)

# Evolve with hydro-gravity bridge
def gravhydro(Mdisk=1e-3, Ndisk=1000):
    rho = float(rho_0(Mdisk*SMBH0.mass, r0).value_in(units.g/units.cm**3))
    print("Running Bridge (Mdisk="+str(round(Mdisk,3))+"Msmbh, rho="+str(rho)+"g/cm3, Ndisk="+str(Ndisk)+")......")
    
    bodies = Particles(0)
    bodies.add_particles(SMBH0)
    bodies.add_particles(BBH0)

    SMBH = bodies[bodies.name=="SMBH"]
    BBH = bodies[bodies.name!="SMBH"]
    BBH = BBH[BBH.name!="gas disk"]

    Nbody = SMBH+BBH

    gravityA = ph4(converter1, number_of_workers=32)
    gravityA.particles.add_particles(Nbody)
    channel = {"from_BHs": bodies.new_channel_to(gravityA.particles),
                "to_BHs": gravityA.particles.new_channel_to(bodies)}
    
    gasdisk, Pinner2, converter2 = make_gasdisk_around_SMBH(SMBH0, Mdisk=Mdisk, Ndisk=Ndisk)

    hydro = Fi(converter2, mode="openmp", workers=32)
    hydro.parameters.use_hydro_flag = True
    hydro.parameters.radiation_flag = False
    hydro.parameters.gamma = 1
    hydro.parameters.isothermal_flag = True
    hydro.parameters.integrate_entropy_flag = False
    hydro.parameters.timestep = 5e2 | units.yr
    hydro.parameters.verbosity = 0
    hydro.parameters.eps_is_h_flag = False    # h_smooth is constant
    eps = 10 | units.au
    hydro.parameters.gas_epsilon = eps
    hydro.parameters.sph_h_const = eps
    
    bodies.add_particles(gasdisk)
    gasdisk = bodies[bodies.name=="gas disk"]

    hydro.particles.add_particles(gasdisk)
    channel.update({"from_gasdisk": gasdisk.new_channel_to(hydro.particles)})
    channel.update({"to_gasdisk": hydro.particles.new_channel_to(gasdisk)})
    
    # Bridge gravity and hydro
    gravhydro = bridge.Bridge(use_threading=False) #, method=SPLIT_4TH_S_M4)
    gravhydro.add_system(gravityA, (hydro,))
    gravhydro.add_system(hydro, (gravityA,))

    # Stopping condition of collision 
    stopping_condition = gravityA.stopping_conditions.collision_detection
    stopping_condition.enable()
    collision_radius_multiplication_factor = 1e5

    Nbody.radius *= collision_radius_multiplication_factor

    # Bridge timestep
    gravhydro.timestep = 5e3 | units.yr

    model_time = 0 | units.Myr
    end_time = 10 | units.Myr
    dt = 1e-2 | units.Myr
    x = 0
    
    BBH_hg_info = []

    gravity_initial_total_energy = gravityA.get_total_energy() + hydro.get_total_energy()

    while model_time < end_time:    

        ProgressBar(end_time.value_in(units.Myr), model_time.value_in(units.Myr))
        gravhydro.evolve_model(model_time)
        # collision detection
        resolve_collision(stopping_condition, gravityA, Nbody) 

        channel["to_BHs"].copy()
        channel["to_gasdisk"].copy()

        if model_time.value_in(units.Myr)>=x:
            x+=0.5

            BBH_hg_info_time = []
            for i in range(N_BBH):
                name1, name2 = BBH[i*2:i*2+2].name
                mass1, mass2 = BBH[i*2:i*2+2].mass.value_in(units.MSun)
                a = BBH[i*2].position - BBH[i*+2+1].position
                hardness = get_hardness(BBH[i*2:i*2+2], Nbody)
                a = relative(BBH[i*2].position, BBH[i*2+1].position).value_in(units.au)
                r = relative(BBH[i*2:i*2+2].center_of_mass(), SMBH.position).value_in(units.pc)
                v1_v2 = relative(BBH[i*2].velocity, BBH[i*2+1].velocity).value_in(units.km/units.s)

                BBH_hg_info_time.append([model_time.value_in(units.Myr), 
                                      name1, name2, mass1, mass2, hardness, a, r, v1_v2])

            BBH_hg_info.append(BBH_hg_info_time)

        model_time += dt

    gravityA.stop()
    hydro.stop()
    
    return np.array(BBH_hg_info)