import sys
from scipy.special import *
import math
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def DiffCoeffWilkeLee(p,T,M_fuel,nnC,nH,nO,nCl,nDB,NBP):

  M_air = 28.96
  M_AB = 2.0/((1/M_air) + (1/M_fuel))
  
  # V_b = liquid molar volume @ normal boiling point 
  # tabled values according to Schroeder (Chaper 4):
  V_b_fuel =  7.0 * ( nC + nH + nO + nDB) + 24.5*nCl

  print('V_b',V_b_fuel)

  sigma_air = 3.617
  sigma_fuel = 1.18*V_b_fuel**(1.0/3.0)  
  sigma_AB = 0.5* (sigma_air + sigma_fuel)

  epsilon_air = 97.0 # air /K
  epsilon_fuel = 1.15*NBP # normal boiling point /K
  epsilon_AB = (epsilon_air * epsilon_fuel)**(0.5) # /K

  T_star = T / epsilon_AB

  A = 1.06036 
  B = 0.15610
  C = 0.19300
  D = 0.47635
  E = 1.03587
  F = 1.52996
  G = 1.76474
  H = 3.89411

  Omega_D = (A/(T_star**B)) + (C/(np.exp(D * T_star))) + (E/(np.exp(F * T_star))) + (G/(np.exp(H * T_star)))

  DiffusionCoefficientWilkeLee = ((3.03 - 0.98/(M_AB**(0.5)))*(1E-3) * T**(1.5))/((p/(100000.0))*(M_AB**(0.5)) * (sigma_AB**2)*Omega_D)/(10000.0)

  return DiffusionCoefficientWilkeLee

#SAE colors
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

dashList = [(5,2),(2,5),(4,10),(3,3,2,2),(5,2,20,2)]


p = 100000
T = 298
M_fuel = 76.5
nC = 3
nH = 5
nO = 0
nCl = 1
nDB = 1
NBP = 318.3

p = 6000000
T = 900
M_fuel = 170.338
nC = 12
nH = 26
nO = 0
nCl = 0
nDB = 0
NBP = 489.47

DAB = DiffCoeffWilkeLee(p,T,M_fuel,nC,nH,nO,nCl,nDB,NBP)
print(DAB*10000)

'''
fsz=24
plt.rc('text', usetex=True)
font = {'family':'serif','size':fsz, 'serif': ['Computer Modern Roman']}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':fsz})


nd = 10000
Red = np.arange(10,nd)

CD_SN = 24.0/Red*(1.0+0.15*Red**0.687)
CD_WA = 24.0/Red*(1.0+1.0/6.0*Red**(2.0/3.0))

CD_SN[Red>1000] = 0.424
CD_WA[Red>1000] = 0.424

f=open('DragModelComparison.dat','w')
np.savetxt(f,np.array([Red,CD_SN,CD_WA]).T, delimiter=" ",fmt='%1.16E')

plt.xlabel(r'Re$_{\mathrm{d}}$')
plt.ylabel(r'$C_{\mathrm{drag}}$')
plt.semilogx(Red, CD_SN,ls='--', lw =2, color=SAEdblue,label='Schiller \& Naumann (1935)')
plt.semilogx(Red, CD_WA,ls='-', lw=2, color=SAEred,label='Wallis (1969)')
plt.legend()
plt.savefig('DragModelComparison.pdf', bbox_inches='tight')

plt.show()

'''
