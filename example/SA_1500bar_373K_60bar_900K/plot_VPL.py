import sys
import numpy as np
from matplotlib import pyplot as plt

SPL1d = np.loadtxt("PL_1d.out", comments="x", unpack=False)
SPL2d = np.loadtxt("PL_2d.out", comments="x", unpack=False)
SPL2dsd = np.loadtxt("PL_dsd_2d.out", comments="x", unpack=False)
SPL4d = np.loadtxt("PL_4d.out", comments="x", unpack=False)
SPL8d = np.loadtxt("PL_8d.out", comments="x", unpack=False)
EVPL = np.loadtxt("SprayA_Evap.dat", comments="x", unpack=False)
plt.plot(SPL1d[:,0],SPL1d[:,2], color='red', lw=2, label='Sim VPL: grid dz=1xD')
plt.plot(SPL2d[:,0],SPL2d[:,2], color='red', ls='--', lw=2, label='Sim VPL: grid dz=2xD')
plt.plot(SPL2dsd[:,0],SPL2dsd[:,2], color='red', ls='--', lw=3, label='Sim VPL: grid dz=2xD dsd')
plt.plot(SPL4d[:,0],SPL4d[:,2], color='red', ls='-.', lw=2, label='Sim VPL: grid dz=4xD')
plt.plot(SPL8d[:,0],SPL8d[:,2], color='red', lw=1, label='Sim VPL: grid dz=8xD')
plt.plot(EVPL[:,0],EVPL[:,1], linestyle="None",marker='o',mfc='None',mec='red',mew=1,markevery=3, label='Exp: VPL')
plt.ylabel('Spray Penetration Length [mm]')
plt.xlabel('time [ms]')
plt.title('Spray Penetration Vs Time')
plt.xlim([0,1.0])
plt.ylim([0,50.0])
legend = plt.legend(loc='lower right', shadow=True)
plt.savefig('VPL_mesh_convergence_dsd.pdf', bbox_inches='tight')
plt.show()
