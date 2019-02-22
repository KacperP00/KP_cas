import sys
import numpy as np
from matplotlib import pyplot as plt

SPL1d = np.loadtxt("PL_1d.out", comments="x", unpack=False)
SPL2d = np.loadtxt("PL_2d.out", comments="x", unpack=False)
SPL2dsd = np.loadtxt("PL_dsd_2d.out", comments="x", unpack=False)
SPL4d = np.loadtxt("PL_4d.out", comments="x", unpack=False)
SPL8d = np.loadtxt("PL_8d.out", comments="x", unpack=False)
ELPL = np.loadtxt("Liq_SprayA_Evap.dat", comments="x", unpack=False)
plt.plot(SPL1d[:,0],SPL1d[:,1], color='blue', lw=2, label='Sim LPL: grid dz=1xD')
plt.plot(SPL2d[:,0],SPL2d[:,1], color='blue', ls='--', lw=2, label='Sim LPL: grid dz=2xD')
plt.plot(SPL2dsd[:,0],SPL2dsd[:,1], color='blue', ls='--', lw=3, label='Sim LPL: grid dz=2xD dsd')
plt.plot(SPL4d[:,0],SPL4d[:,1], color='blue', ls='-.', lw=2, label='Sim LPL: grid dz=4xD')
plt.plot(SPL8d[:,0],SPL8d[:,1], color='blue',  lw=1, label='Sim LPL: grid dz=8xD')
plt.plot(ELPL[:,0],ELPL[:,1], linestyle="None",marker='s',mfc='None',mec='blue',mew=1,markevery=3, label='Exp: LPL')
plt.ylabel('Spray Penetration Length [mm]')
plt.xlabel('time [ms]')
plt.title('Spray Penetration Vs Time')
plt.xlim([0,1.0])
plt.ylim([0,20.0])
legend = plt.legend(loc='lower right', shadow=True)
plt.savefig('LPL_mesh_convergence_dsd.pdf', bbox_inches='tight')
plt.show()
