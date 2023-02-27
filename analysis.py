import numpy as np
import matplotlib.pyplot as plt
from amuse.support import io
from amuse.units import units

BBH0 = io.read_set_from_file(
       'BBH0.txt', 'txt',
       attribute_types = (units.MSun, units.pc, units.kms, units.kms, units.kms, units.pc, units.pc, units.pc),
       attribute_names= ("mass", "radius", "vx", "vy", "vz", "x", "y", "z")
)
N_BBH = int(len(BBH0)/2)

BBH_bridge = np.load("Results/Bridge/BBH_bridge.npy")
BBH_bridge_m01 = np.load("Results/Bridge/BBH_bridge_m01.npy")
BBH_bridge_m10 = np.load("Results/Bridge/BBH_bridge_m10.npy")
BBH_nbody = np.load("Results/Nbody/BBH_nbody.npy")
BBH_nbody_m01 = np.load("Results/Nbody/BBH_nbody_m01.npy")
BBH_nbody_m10 = np.load("Results/Nbody/BBH_nbody_m10.npy")
BBH_ref_info = np.load("Results/Reference/BBH_ref_info.npy")


harder_BBH_bridge = 0
softer_BBH_bridge = 0
not_BBH_bridge = 0
much_harder_BBH_bridge = 0

harder_BBH_bridge_m01 = 0
softer_BBH_bridge_m01 = 0
not_BBH_bridge_m01 = 0
much_harder_BBH_bridge_m01 = 0

harder_BBH_bridge_m10 = 0
softer_BBH_bridge_m10 = 0
not_BBH_bridge_m10 = 0
much_harder_BBH_bridge_m10 = 0

harder_BBH_nbody = 0
softer_BBH_nbody = 0
not_BBH_nbody = 0
much_harder_BBH_nbody = 0

harder_BBH_nbody_m01 = 0
softer_BBH_nbody_m01 = 0
not_BBH_nbody_m01 = 0
much_harder_BBH_nbody_m01 = 0

harder_BBH_nbody_m10 = 0
softer_BBH_nbody_m10 = 0
not_BBH_nbody_m10 = 0
much_harder_BBH_nbody_m10 = 0

softerness = []
harderness = []
softerness_m01 = []
harderness_m01 = []
softerness_m10 = []
harderness_m10 = []
softerness_n = []
harderness_n = []
softerness_n_m01 = []
harderness_n_m01 = []
softerness_n_m10 = []
harderness_n_m10 = []


# Bridge_Mdisk=Msmbh
for i in range(10):
    for j in range(N_BBH):
        hardness = np.array(BBH_bridge)[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_bridge += 1
        elif d_hardness[-1] < 0:
            softer_BBH_bridge += 1
            softerness.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_bridge += 1  
            harderness.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_bridge += 1  

# Bridge_Mdisk=0.1Msmbh
for i in range(10):
    for j in range(N_BBH):
        hardness = np.array(BBH_bridge_m01)[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_bridge_m01 += 1
        elif d_hardness[-1] < 0:
            softer_BBH_bridge_m01 += 1
            softerness_m01.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_bridge_m01 += 1  
            harderness_m01.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_bridge_m01 += 1 

# Bridge_Mdisk=10Msmbh                
for i in range(10):
    for j in range(N_BBH):
        hardness = np.array(BBH_bridge_m10)[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_bridge_m10 += 1
        elif d_hardness[-1] < 0:
            softer_BBH_bridge_m10 += 1
            softerness_m10.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_bridge_m10 += 1  
            harderness_m10.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_bridge_m10 += 1 

# Nbody_Mmin=Msun                
for i in range(10):
    for j in range(N_BBH):
        hardness = np.array(BBH_nbody)[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_nbody += 1
        elif d_hardness[-1] < 0:
            softer_BBH_nbody += 1
            softerness_n.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_nbody += 1  
            harderness_n.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_nbody += 1 

# Nbody_Mmin=0.1Msun                
for i in range(10):
    for j in range(N_BBH):
        hardness = np.array(BBH_nbody_m01)[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_nbody_m01 += 1
        elif d_hardness[-1] < 0:
            softer_BBH_nbody_m01 += 1
            softerness_n_m01.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_nbody_m01 += 1  
            harderness_n_m01.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_nbody_m01 += 1 
                
# Nbody_Mmin=10Msun                
for i in range(10):
    for j in range(N_BBH):
        hardness = np.array(BBH_nbody_m10)[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_nbody_m10 += 1
        elif d_hardness[-1] < 0:
            softer_BBH_nbody_m10 += 1
            softerness_n_m10.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_nbody_m10 += 1  
            harderness_n_m10.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_nbody_m10 += 1 
                

plt.figure()
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(softerness_n,axis=0).astype('float'), linestyle="dashed", c="r")
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(harderness_n,axis=0).astype('float'), c="r", label=r"$M_{min}=M_{sun}$")
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(softerness_n_m01,axis=0).astype('float'), linestyle="dashed", c="g")
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(harderness_n_m01,axis=0).astype('float'), c="g", label=r"$M_{min}=0.1M_{sun}$")
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(softerness_n_m10,axis=0).astype('float'), linestyle="dashed", c="b")
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(harderness_n_m10,axis=0).astype('float'), c="b", label=r"$M_{min}=10M_{sun}$")
plt.plot([-1,10],[0,0], linestyle="dotted", c="gray")
plt.annotate("Harder", xy=(0.01,0.95), xycoords='axes fraction')
plt.annotate("Softer", xy=(0.01,0.01), xycoords='axes fraction')
plt.xlabel('Time [Myr]', size=12)
plt.ylabel(r'$\Delta$Hardness', size=12)
plt.legend(bbox_to_anchor=(1.0,1.0))
plt.xlim([0,9.5])
plt.title("Average hardness variation of BBHs in gas-free galactic nuclei")
plt.savefig("Figures/hardness_nbody.png", dpi=200, bbox_inches='tight')
plt.show()



plt.figure()
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(softerness,axis=0).astype('float'), c="r", linestyle="dashed")
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(softerness_m01,axis=0).astype('float'), c="g", linestyle="dashed")
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(softerness_m10,axis=0).astype('float'), c="b", linestyle="dashed")


plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(harderness,axis=0).astype('float'), c="r", label=r"M$_{disk}$=M$_{SMBH}$")
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(harderness_m01,axis=0).astype('float'), c="g", label=r"M$_{disk}$=0.1M$_{SMBH}$")
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(harderness_m10,axis=0).astype('float'), c="b", label=r"M$_{disk}$=10M$_{SMBH}$")


plt.plot([-1,10],[0,0], linestyle="dotted", c="gray")
plt.annotate("Harder", xy=(0.01,0.95), xycoords='axes fraction')
plt.annotate("Softer", xy=(0.01,0.01), xycoords='axes fraction')
plt.xlabel('Time [Myr]', size=12)
plt.ylabel(r'$\Delta$Hardness', size=12)
plt.legend(bbox_to_anchor=(1.0,1.0))
plt.xlim([0,9.5])
plt.title("Average hardness variation of BBHs in gaseous galactic nuclei")
plt.savefig("Figures/hardness_bridge.png", dpi=200, bbox_inches='tight')
plt.show()




r5 = BBH_ref_info[0,:,7].astype("float")<5
r10 = BBH_ref_info[0,:,7].astype("float")>5

harder_BBH_nbody = 0
softer_BBH_nbody = 0
not_BBH_nbody = 0
much_harder_BBH_nbody = 0

harder_BBH_nbody_r5 = 0
softer_BBH_nbody_r5 = 0
not_BBH_nbody_r5 = 0
much_harder_BBH_nbody_r5 = 0

harder_BBH_nbody_r10 = 0
softer_BBH_nbody_r10 = 0
not_BBH_nbody_r10 = 0
much_harder_BBH_nbody_r10 = 0

softerness = []
harderness = []
softerness_r5 = []
harderness_r5 = []
softerness_r10 = []
harderness_r10 = []

BBH_nbody_r5 = np.array(BBH_nbody)[:,:,r5,:]
BBH_nbody_r10 = np.array(BBH_nbody)[:,:,r10,:]
BBH_ref_info_r5 = BBH_ref_info[:,r5,:]
BBH_ref_info_r10 = BBH_ref_info[:,r10,:]


for i in range(10):
    for j in range(N_BBH):
        hardness = np.array(BBH_nbody)[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_nbody += 1
        elif d_hardness[-1] < 0:
            softer_BBH_nbody += 1
            softerness.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_nbody += 1  
            harderness.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_nbody += 1  
                
for i in range(10):
    for j in range(BBH_nbody_r5.shape[2]):
        hardness = BBH_nbody_r5[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info_r5[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_nbody_r5 += 1
        elif d_hardness[-1] < 0:
            softer_BBH_nbody_r5 += 1
            softerness_r5.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_nbody_r5 += 1  
            harderness_r5.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_nbody_r5 += 1 
                
for i in range(10):
    for j in range(BBH_nbody_r10.shape[2]):
        hardness = BBH_nbody_r10[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info_r10[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_nbody_r10 += 1
        elif d_hardness[-1] < 0:
            softer_BBH_nbody_r10 += 1
            softerness_r10.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_nbody_r10 += 1  
            harderness_r10.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_nbody_r10 += 1 
                

softerness = np.array(softerness)
harderness = np.array(harderness)


plt.figure()
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(softerness,axis=0).astype('float'), c="r", linestyle="dashed")
plt.plot(BBH_ref_info_r5[:,0,0].astype('float'), np.mean(softerness_r5,axis=0).astype('float'), c="g", linestyle="dashed")
plt.plot(BBH_ref_info_r10[:,0,0].astype('float'), np.mean(softerness_r10,axis=0).astype('float'), c="b", linestyle="dashed")
 
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(harderness,axis=0).astype('float'), c="r", label="BBHs")
plt.plot(BBH_ref_info_r5[:,0,0].astype('float'), np.mean(harderness_r5,axis=0).astype('float'), c="g", label="BBHs(r < 5 pc)")
plt.plot(BBH_ref_info_r10[:,0,0].astype('float'), np.mean(harderness_r10,axis=0).astype('float'), c="b", label="BBHs(5 pc < r < 10 pc)")

plt.plot([-1,10],[0,0], linestyle="dotted", c="gray")
plt.annotate("Harder", xy=(0.01,0.95), xycoords='axes fraction')
plt.annotate("Softer", xy=(0.01,0.01), xycoords='axes fraction')
plt.xlabel('Time [Myr]', size=12)
plt.ylabel(r'$\Delta$Hardness', size=12)
plt.legend(bbox_to_anchor=(1.0,1.0))
plt.xlim([0,9.5])
plt.title("Average hardness variation of BBHs in gas-free galactic nuclei")
plt.savefig("Figures/hardness_nbody_m=1.png", dpi=200, bbox_inches='tight')
plt.show()


harder_BBH_bridge = 0
softer_BBH_bridge = 0
not_BBH_bridge = 0
much_harder_BBH_bridge = 0

harder_BBH_bridge_r5 = 0
softer_BBH_bridge_r5 = 0
not_BBH_bridge_r5 = 0
much_harder_BBH_bridge_r5 = 0

harder_BBH_bridge_r10 = 0
softer_BBH_bridge_r10 = 0
not_BBH_bridge_r10 = 0
much_harder_BBH_bridge_r10 = 0

softerness = []
harderness = []
softerness_r5 = []
harderness_r5 = []
softerness_r10 = []
harderness_r10 = []

BBH_bridge_r5 = np.array(BBH_bridge)[:,:,r5,:]
BBH_bridge_r10 = np.array(BBH_bridge)[:,:,r10,:]
BBH_ref_info_r5 = BBH_ref_info[:,r5,:]
BBH_ref_info_r10 = BBH_ref_info[:,r10,:]


for i in range(10):
    for j in range(N_BBH):
        hardness = np.array(BBH_bridge)[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_bridge += 1
        elif d_hardness[-1] < 0:
            softer_BBH_bridge += 1
            softerness.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_bridge += 1  
            harderness.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_bridge += 1  
                
for i in range(10):
    for j in range(BBH_bridge_r5.shape[2]):
        hardness = BBH_bridge_r5[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info_r5[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_bridge_r5 += 1
        elif d_hardness[-1] < 0:
            softer_BBH_bridge_r5 += 1
            softerness_r5.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_bridge_r5 += 1  
            harderness_r5.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_bridge_r5 += 1 
                
for i in range(10):
    for j in range(BBH_bridge_r10.shape[2]):
        hardness = BBH_bridge_r10[i,:,j,5].astype('float')
        d_hardness = hardness - BBH_ref_info_r10[:,j,5].astype('float')
        if hardness[-1] < 0:
            not_BBH_bridge_r10 += 1
        elif d_hardness[-1] < 0:
            softer_BBH_bridge_r10 += 1
            softerness_r10.append(d_hardness)
        elif d_hardness[-1] > 0:
            harder_BBH_bridge_r10 += 1  
            harderness_r10.append(d_hardness)
            if hardness[-1] > 10:
                much_harder_BBH_bridge_r10 += 1 
                

softerness = np.array(softerness)
harderness = np.array(harderness)


plt.figure()
plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(softerness,axis=0).astype('float'), c="r", linestyle="dashed")
plt.plot(BBH_ref_info_r5[:,0,0].astype('float'), np.mean(softerness_r5,axis=0).astype('float'), c="g", linestyle="dashed")
plt.plot(BBH_ref_info_r10[:,0,0].astype('float'), np.mean(softerness_r10,axis=0).astype('float'), c="b", linestyle="dashed")

plt.plot(BBH_ref_info[:,0,0].astype('float'), np.mean(harderness,axis=0).astype('float'), c="r" , label="BBHs")
plt.plot(BBH_ref_info_r5[:,0,0].astype('float'), np.mean(harderness_r5,axis=0).astype('float'), c="g", label="BBHs(r < 5 pc)")
plt.plot(BBH_ref_info_r10[:,0,0].astype('float'), np.mean(harderness_r10,axis=0).astype('float'), c="b", label="BBHs(5 pc < r < 10 pc)")

plt.plot([-1,10],[0,0], linestyle="dotted", c="gray")
plt.annotate("Harder", xy=(0.01,0.95), xycoords='axes fraction')
plt.annotate("Softer", xy=(0.01,0.01), xycoords='axes fraction')
plt.xlabel('Time [Myr]', size=12)
plt.ylabel(r'$\Delta$Hardness', size=12)
plt.legend(bbox_to_anchor=(1.0,1.0))
plt.xlim([0,9.5])
plt.title("Average hardness variation of BBHs in gaseous galactic nuclei")
plt.savefig("Figures/hardness_bridge_m=1.png", dpi=200, bbox_inches='tight')
plt.show()