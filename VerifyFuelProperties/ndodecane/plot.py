import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import rc
import math
from scipy.special import *

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

data0 = np.loadtxt('Patm_Tatm/LiquidFuelProperties.dat',comments="#", unpack=False, skiprows=2)
data1 = np.loadtxt('Pcrit_Tcrit/LiquidFuelProperties.dat',comments="#", unpack=False, skiprows=2)
data2 = np.loadtxt('constant_Lv/LiquidFuelProperties.dat',comments="#", unpack=False, skiprows=2)
data3 = np.loadtxt('Pcrit_const_Lv/LiquidFuelProperties.dat',comments="#", unpack=False, skiprows=2)

#rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
fs = 16
params = {'legend.fontsize': fs,
         'axes.labelsize': fs,
         'axes.titlesize': fs,
         'xtick.labelsize': fs,
         'ytick.labelsize': fs}
mpl.rcParams.update(params)
plt.grid(b=True, which='major', linestyle=':')
#plt.yscale('log',basey=10)
#plt.xscale('log',basey=10)
plt.xlim([300,800])
plt.ylim([0,20])

plt.plot(data0[:,0],data0[:,7]/1e5, lw=2, color=SAEred, label='Daubert Danner')
plt.plot(data0[:,0],data0[:,9]/1e5, lw=2, color=SAEdgreen, label=r'CC: P$_{atm}$')
plt.plot(data1[:,0],data1[:,9]/1e5, lw=2, color=SAEdblue, label=r'CC: P$_{crit}$')
plt.plot(data2[:,0],data2[:,9]/1e5, lw=2, color=SAEblue, label=r'CC: const. LHV')
plt.plot(data3[:,0],data3[:,9]/1e5, lw=2, color=SAEdgray, label=r'CC: P$_{crit}$ const. LHV')

plt.ylabel('Vapor Pressure / bar')
plt.xlabel('Temperature / K')
legend = plt.legend(loc='upper left', shadow=True)
plt.savefig('vapor_pressure.pdf', bbox_inches='tight')
plt.show()

