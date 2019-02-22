import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import rc

#data = np.loadtxt(sys.argv[1],comments="#", delimiter="   ", unpack=False)
data = np.loadtxt(sys.argv[1],comments="#", unpack=False)
DmD2 = np.loadtxt("DmD2_SprayA_Evap.dat",comments="#", unpack=False)


SAEyellow = np.array(np.asarray([255.0,178.0,  1.0])/255.0).tolist()
SAEorange = np.array(np.asarray([234.0,113.0, 37.0])/255.0).tolist()
SAEred    = np.array(np.asarray([220.0, 41.0, 30.0])/255.0).tolist()
SAEblue   = np.array(np.asarray([  1.0,160.0,233.0])/255.0).tolist()
SAEdblue  = np.array(np.asarray([  0.0, 81.0,149.0])/255.0).tolist()
SAEgreen  = np.array(np.asarray([ 46.0,177.0, 53.0])/255.0).tolist()
SAEdgreen = np.array(np.asarray([  0.0,119.0, 61.0])/255.0).tolist()
SAElgray  = np.array(np.asarray([202.0,202.0,200.0])/255.0).tolist()
SAEmgray  = np.array(np.asarray([154.0,155.0,157.0])/255.0).tolist()
SAEdgray  = np.array(np.asarray([ 97.0, 98.0,101.0])/255.0).tolist()
SAEblack  = np.array(np.asarray([  0.0,  0.0,  0.0])/255.0).tolist()

colr = [SAEred, SAEdblue, SAEdgreen, SAEorange, SAEblue, SAEgreen, SAEyellow, SAEblack, SAEdgray]
nclr = max(np.shape(colr))

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
fs = 16
params = {'legend.fontsize': fs,
         'axes.labelsize': fs,
         'axes.titlesize': fs,
         'xtick.labelsize': fs,
         'ytick.labelsize': fs}
mpl.rcParams.update(params)
plt.grid(b=True, which='major', linestyle=':')
plt.ylabel(r'TKE / m$^2$/s$^2$ OR $\chi_g$ / s^{-1}')
plt.xlabel(r'z / mm')
#plt.yscale('log',basey=10)
#plt.xscale('log',basey=10)
#plt.ylim([1,10000])

plt.plot(data[:,0]*90e-3,data[:,14]*0.5*5.73446E+002**2, lw=2, color='blue', label=r'$k_g$')
#plt.plot(data[:,0],data[:,7], lw=2, label='u_g',color='green')
#plt.plot(data[:,0],data[:,16], lw=2, color='blue', label='mu_t_g')
plt.plot(data[:,0]*90e-3,2.*data[:,15]/data[:,14]*data[:,17]*5.73446E+002/90e-6, lw=2, color='red', label=r'$\chi_g$')
#plt.plot(data[:,0]*90e-3,(data[:,0]*90)**(-4),color='black',label=r'$z^{-4}$')

legend = plt.legend(loc='upper right', shadow=True)

plt.twinx().plot(data[:,0]*90e-3,data[:,15]*0.5*5.73446E+002**3/90e-6, lw=2, color='green', label=r'$\varepsilon_g$')

plt.ylabel(r'$\varepsilon_g$ / m$^2$/s$^3$')
#plt.xlabel('z/D [-]')
#plt.xlim([0,40])
#plt.ylim([0,5e8])
legend = plt.legend(loc='lower right', shadow=True)
plt.savefig('turbulence.pdf', bbox_inches='tight')
plt.show()


plt.clf()
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
fs = 16
params = {'legend.fontsize': fs,
         'axes.labelsize': fs,
         'axes.titlesize': fs,
         'xtick.labelsize': fs,
         'ytick.labelsize': fs}
mpl.rcParams.update(params)
plt.grid(b=True, which='major', linestyle=':')
plt.ylabel(r'Z$_{mix}$ / -')
plt.xlabel(r'z / mm')
#plt.ylim([0,0.25])

plt.plot(data[:,0]*90e-3,(3.475402137*data[:,3]-0.232*data[:,4]+0.232)/(3.475402137+0.232), lw=2, color='blue', label=r'Z$_{mix}$')
legend = plt.legend(loc='upper left', shadow=True)

plt.twinx().plot(data[:,0]*90e-3,data[:,17], lw=2, color='red', label=r'Z$_{var}$')

plt.ylabel(r'Z$_{var}$ / -')
#plt.xlabel('z/D [-]')
#plt.xlim([0,40])
#plt.ylim([0,3.0e-5])
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('mixture_frac_var.pdf', bbox_inches='tight')
plt.show()


plt.clf()
plt.plot(data[:,0],data[:,2], lw=2, label='Y_l')
plt.plot(data[:,0],data[:,3], lw=2, label='Y_v')
plt.plot(data[:,0],data[:,4], lw=2, label='Y_a')
#plt.plot(data[:,0],data[:,14], ls='None', marker='.', lw=2, label='Y_ref')
plt.ylabel('Mass Fraction')
plt.xlabel('z/D [-]')
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('mass_fraction.pdf', bbox_inches='tight')
plt.show()


'''
plt.clf()
plt.plot(data[:,0],data[:,14], lw=2, label='k_g')
plt.plot(data[:,0],data[:,15], ls='None', marker='.', lw=2, label=r'\epsilon_g')
#plt.plot(data[:,0],data[:,16], lw=2, label='mu_t_g')
#plt.plot(data[:,0],0.09*data[:,1]*(data[:,4]+data[:,3])*data[:,14]**2/data[:,15], lw=2, label='mu_t_g_m')
plt.ylabel('Turbulence')
plt.xlabel('z/D [-]')
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('turb_keps.pdf', bbox_inches='tight')
plt.show()



plt.clf()
#plt.hist(data[:100,14])
#plt.plot(data[:,0],data[:,1]*(data[:,2]+data[:,3]+data[:,4])*data[:,13]**2)
plt.plot(data[:,0],data[:,15],label='1')
plt.plot(data[:,0],data[:,16],label='2')
plt.plot(data[:,0],data[:,18],label='4')
legend = plt.legend(loc='upper left', shadow=True)
plt.show()

'''
plt.clf()
plt.plot(data[:,0],data[:,1], lw=2, label='rho')
plt.plot(data[:,0],data[:,6], lw=2, label='u_l')
plt.plot(data[:,0],data[:,7], lw=2, label='u_g')
plt.plot(data[:,0],data[:,11], lw=2, label='T_d')
plt.plot(data[:,0],data[:,12], lw=2, label='T_g')
plt.ylabel('Flow Variables')
plt.xlabel('z/D [-]')
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('flow_variables.pdf', bbox_inches='tight')
plt.show()
'''

plt.clf()
plt.plot(data[:,0],data[:,3]**2*(0.039+24.8*(63889.7*data[:,6]*2*data[:,13]/31.0)**(-0.647)), lw=2, label='b')
#plt.plot(data[:,0],(63889.7*data[:,6]*data[:,13]/31.0), lw=2, label='Reg_z')
plt.ylabel('Spray boundary')
plt.xlabel('z/D [-]')
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('spray_bound.pdf', bbox_inches='tight')
plt.show()

'''
plt.clf()
ax = plt.gca()
ax2 = ax.twinx()
plt.axis('normal')
ax.plot(data[:,0],data[:,8], color='red', lw=2, label='<dm>')
ax.plot(DmD2[:,0],DmD2[:,1],linestyle="None",marker='s',mfc='None',mec='red',mew=1,markevery=1, label='Dm Exp')
ax.set_ylabel("Mean of Diameter",fontsize=14,color='red')
ax2.plot(data[:,0],data[:,9], color='blue', lw=2, label='<d2>')
ax2.set_ylabel("Second Moment of Diameter",fontsize=14,color='blue')
ax2.plot(DmD2[:,0],DmD2[:,2],linestyle="None",marker='s',mfc='None',mec='blue',mew=1,markevery=1, label='D2 Exp')
#ax.set_ylim(ymax=100)
#ax.set_xlim(xmax=100)
#ax.grid(True)
plt.title("Droplet Diameter", fontsize=20,color='black')
ax.set_xlabel('z/D [-]', fontsize=14, color='b')
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('droplet_diameter.pdf', bbox_inches='tight')
plt.show()

plt.clf()
SPL = np.loadtxt("PL_dsd_2d.out", comments="#", unpack=False)
ELPL = np.loadtxt("Liq_SprayA_Evap.dat", comments="#", unpack=False)
EVPL = np.loadtxt("Vap_SprayA_Evap.dat", comments="#", unpack=False)
#CVG = np.loadtxt("spray_ecn_all.out", comments="#", unpack=False)
plt.plot(SPL[:,0],SPL[:,1], color='blue', lw=2, label='Sim: LPL')
#plt.plot(CVG[:,0]*1000,CVG[:,6]*1000,linestyle="None",marker='s',mfc='None',mec='blue',mew=1,markevery=5, label='CVG: LPL')
plt.plot(ELPL[:,0],ELPL[:,1], linestyle="None",marker='s',mfc='None',mec='blue',mew=1,markevery=5, label='Exp: LPL')
#plt.errorbar(ELPL[:,0],ELPL[:,1],yerr=ELPL[:,3], linestyle="None",marker='s',mfc='None',mec='blue',mew=1,markevery=5, label='Exp: LPL')
#plt.plot(ELPL[:,0],ELPL[:,1], linestyle="-",lw=2,color='blue', label='Exp: LPL')
plt.fill_between(ELPL[:,0],ELPL[:,1]-ELPL[:,3],ELPL[:,1]+ELPL[:,3],alpha=0.2, color='blue', antialiased=True)
plt.plot(SPL[:,0],SPL[:,2], color='red', lw=2, label='Sim: VPL')
#plt.plot(CVG[:,0]*1000,CVG[:,7]*1000,linestyle="None",marker='s',mfc='None',mec='red',mew=1,markevery=5, label='CVG: VPL')
plt.plot(EVPL[:,0],EVPL[:,1], linestyle="None",marker='o',mfc='None',mec='red',mew=1,markevery=5, label='Exp: VPL')
#plt.errorbar(EVPL[:,0],EVPL[:,1],yerr=EVPL[:,3], linestyle="None",marker='o',mfc='None',mec='red',mew=1,markevery=5, label='Exp: VPL')
#plt.plot(EVPL[:,0],EVPL[:,1], linestyle="-",lw=2,color='red',label='Exp: VPL')
plt.fill_between(EVPL[:,0],EVPL[:,1]-EVPL[:,3],EVPL[:,1]+EVPL[:,3],alpha=0.2, color='red', antialiased=True)
plt.ylabel('Spray Penetration Length [mm]')
plt.xlabel('time [ms]')
plt.title('Spray Penetration Vs Time')
plt.xlim([0,1.0])
plt.ylim([0,50.0])
legend = plt.legend(loc='upper left', shadow=True)
plt.savefig('PL_sim_exp.pdf', bbox_inches='tight')
plt.show()

'''
plt.clf()
SPL1d = np.loadtxt("PL_1d.out", comments="#", unpack=False)
SPL2d = np.loadtxt("PL_2d.out", comments="#", unpack=False)
SPL4d = np.loadtxt("PL_4d.out", comments="#", unpack=False)
SPL8d = np.loadtxt("PL_8d.out", comments="#", unpack=False)
ELPL = np.loadtxt("Liq_SprayA_Evap.dat", comments="#", unpack=False)
plt.plot(SPL1d[:,0],SPL1d[:,1], color='blue', lw=2, label='Sim LPL: grid dz=1xD')
plt.plot(SPL2d[:,0],SPL2d[:,1], color='blue', ls='--', lw=2, label='Sim LPL: grid dz=2xD')
plt.plot(SPL4d[:,0],SPL4d[:,1], color='blue', ls='-.', lw=2, label='Sim LPL: grid dz=4xD')
plt.plot(SPL8d[:,0],SPL8d[:,1], color='blue',  lw=1, label='Sim LPL: grid dz=8xD')
plt.plot(ELPL[:,0],ELPL[:,1], linestyle="None",marker='s',mfc='None',mec='blue',mew=1,markevery=3, label='Exp: LPL')
plt.ylabel('Spray Penetration Length [mm]')
plt.xlabel('time [ms]')
plt.title('Spray Penetration Vs Time')
plt.xlim([0,1.0])
plt.ylim([0,20.0])
legend = plt.legend(loc='lower right', shadow=True)
plt.savefig('LPL_mesh_convergence.pdf', bbox_inches='tight')
plt.show()

plt.clf()
SPL1d = np.loadtxt("PL_1d.out", comments="#", unpack=False)
SPL2d = np.loadtxt("PL_2d.out", comments="#", unpack=False)
SPL4d = np.loadtxt("PL_4d.out", comments="#", unpack=False)
SPL8d = np.loadtxt("PL_8d.out", comments="#", unpack=False)
EVPL = np.loadtxt("Vap_SprayA_Evap.dat", comments="#", unpack=False)
plt.plot(SPL1d[:,0],SPL1d[:,2], color='red', lw=2, label='Sim VPL: grid dz=1xD')
plt.plot(SPL2d[:,0],SPL2d[:,2], color='red', ls='--', lw=2, label='Sim VPL: grid dz=2xD')
plt.plot(SPL4d[:,0],SPL4d[:,2], color='red', ls='-.', lw=2, label='Sim VPL: grid dz=4xD')
plt.plot(SPL8d[:,0],SPL8d[:,2], color='red', lw=1, label='Sim VPL: grid dz=8xD')
plt.plot(EVPL[:,0],EVPL[:,1], linestyle="None",marker='o',mfc='None',mec='red',mew=1,markevery=3, label='Exp: VPL')
plt.ylabel('Spray Penetration Length [mm]')
plt.xlabel('time [ms]')
plt.title('Spray Penetration Vs Time')
plt.xlim([0,1.0])
plt.ylim([0,50.0])
legend = plt.legend(loc='lower right', shadow=True)
plt.savefig('VPL_mesh_convergence.pdf', bbox_inches='tight')
plt.show()
'''
