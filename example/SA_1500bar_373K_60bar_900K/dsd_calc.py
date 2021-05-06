import numpy as np
import matplotlib.pyplot as plt
import csv
import math
from scipy.special import *
'''
## Kelvin-Helmoltz WAVE model DSD calculation
# Fuel properties (n-dodecane)
sigma  = 1.92e-2;                               # Surface tension
rho_l  = 698.3239;                              # Density
visc_l = 5.6946e-4;                             # Dynamic viscosity liquid phase
T_l    = 363;                                   # Fuel Temperature
MW_f = 170.3380;                                # Molecular weight

# Ambient properties
MW_a = 28e-3;                                   # Molecular weight of ambient gas
R = 8.3144598;                                  # Universal gas constant
visc_g = 2.3928e-5;                             # Dynamic viscosity of gas phase
rho_g = 22.8;                                   # Density
T_amb = 900;                                    # Ambient temperature
P_amb = R*rho_g*T_amb/MW_a;                   # Ambient pressure

# Injection parameters
P_inj = 150e6;                                  # Injection pressure
C_d   = 0.89;                                   # Discharge coefficient
u_inj = C_d*np.sqrt(2*(P_inj-P_amb)/rho_l);        # Injection velocity at nozzle exit

Re = 6.46236E+004;
We = 1.12518E+006;
DR = 3.08530E+001;
VR = 1.47592E+001;
B0 = 0.61; 
B1 = 10;
Weg = We/DR/2;
Wel = We/2;
Rel = Re/2;
Z = np.sqrt(Wel)/(Rel);
T = Z*np.sqrt(Weg);
Dn = 90e-6;

m_blob = rho_l*np.pi/6*Dn**3;
M_blob = m_blob*u_inj;
E_blob = 0.5*M_blob*u_inj;

print 'Total Energy: ', E_blob
m_blob0 = m_blob

TotE = 0
ueq = u_inj;
deq = Dn;
i = 0;
di = [];
ui = [];
md = [];
while m_blob>0:
    LambdaKH = (9.02/2.0)*(1.0+0.45*np.sqrt(Z))*(1.0+0.4*T**0.7)/(1.0+0.865*Weg**1.67)**0.6*deq;
    di.append(B0*LambdaKH);
    ui.append(ueq);
    m_di = rho_l*np.pi/6.0*di[i]**3;
    md.append(m_di/m_blob0);
    m_blob = m_blob-m_di;
    if m_blob<=0:
        break;
    M_blob = M_blob-m_di*ueq;
    deq = (6.0*m_blob/np.pi/rho_l)**(1.0/3.0);
    ueq = M_blob/m_blob;
    E_blob = 0.5*M_blob*ueq;
    TotE = TotE + 0.5*m_di*ueq**2;
    print di[i], ui[i],md[i]
    dratio = deq/Dn;
    uratio = ueq/u_inj;
    uratio2 = uratio**2;
    #print dratio,uratio,uratio2
    Weg = Weg*dratio*uratio2;
    Wel = Wel*dratio*uratio2;
    Rel = Rel*dratio*uratio;
    Z = np.sqrt(Wel)/(Rel);
    T = Z*np.sqrt(Weg);
    i = i + 1;
'''
'''
with open('dsd.csv', 'w') as f:
    i = 0
    for val in di:
        f.write("%1.16E,%1.16E\n" % val,md[i])
	i = i + 1
'''
'''
hist, bin_edges = np.histogram(np.array(di)/Dn,bins=100000,density=True)
 
#plt.hist(np.array(md), bins=[0,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1])
#plt.xlim([0,4])
plt.plot(bin_edges[:-1],hist)
plt.show()
'''

def bisection(ratio):
        tol = 1e-6
        ql0 = 1.0
        qh0 = 10.0
        iter = 0
        err = 1.0
        ql = ql0
        qh = qh0
        qs = 0.5*(ql+qh)
        while ( err > tol):
                iter = iter + 1
                fs = gamma(1.0+2.0/qs)/(gamma(1.0+1.0/qs))**2 - ratio
                fl = gamma(1.0+2.0/ql)/(gamma(1.0+1.0/ql))**2 - ratio
                fh = gamma(1.0+2.0/qh)/(gamma(1.0+1.0/qh))**2 - ratio
                if (fl*fs < 0.0):
                        qh = qs; qs = 0.5*(ql+qh)
                elif (fs*fh < 0.0):
                        ql = qs; qs = 0.5*(ql+qh)
                err = abs(fs)
        q_final = qs
        return q_final


def computeDSDrr(Dm,D2):
        L = 10*Dm
        nd = 200
        h = L/nd
        D = h*np.arange(0,nd)
        r_min = 1.05
        r_max = 1.9
        dsd = np.zeros(nd)
        if (Dm > 0.0 and D2 > 0.0):
                ratio = D2/Dm**2
                if (ratio > r_min):
                        ratio = min(r_max,ratio)
                        q = bisection(ratio)
                else:
                        q = 5.1334
                Dbar = Dm/gamma(1.0+1.0/q)
                f = (q/Dbar**q)*D**(q-1.0)*np.exp(-(D/Dbar)**q)
                norm = sum(f*h)
                if (norm > 0.0) :
                        dsd[:] = f/norm
        return [D, dsd]

def computeDSDrl(Dm,D2):
        L = 10*Dm
        nd = 200
        h = L/nd
        D = h*np.arange(0,nd)
        r_min = 1.05
        r_max = 1.9
        dsd = np.zeros(nd)
        if (Dm > 0.0 and D2 > 0.0):
                ratio = D2/Dm**2
                if (ratio > r_min):
                        ratio = min(r_max,ratio)
                        q = bisection(ratio)
                else:
                        q = 5.1334
                Dbar = Dm/gamma(1.0+1.0/q)
                f = (q/Dbar**q)*D**(q-1.0)*np.exp(-(D/Dbar)**q)
                norm = sum(f*h)
                if (norm > 0.0) :
                        dsd[:] = f/norm
        D[0] = 1e-10
        mu = np.log(Dm**2/np.sqrt(D2))
        sigma2 = np.log(D2/Dm**2)
        pi = math.acos(-1.0)
        B = 1.0/np.sqrt(2.0*pi*sigma2)
        dsdln = np.zeros(nd)
        if (Dm > 0.0 and D2 > 0.0):
                f = B*np.exp(-((np.log(D)-mu)**2/2/sigma2))/D
                norm = sum(f*h)
                if (norm > 0.0) :
                        dsdln[:] = f/norm
        peak = max(dsdln)
        location = np.argmax(dsdln,0)
        #print peak,location
        weight = dsdln/peak
        for i in np.arange(1,max(np.shape(dsdln))):
                if (weight[i-1] == 1.0 ):
                        weight[i] = 1.0
        f = (1.0-weight)*dsdln + weight*dsd
        norm = sum(f*h)
        if (norm > 0.0) :
                dsd[:] = f/norm
        return [D, dsd]

def computeDSDgauss(Dm,D2):
        L = 10*Dm
        nd = 200
        h = L/nd
        D = h*np.arange(0,nd)
        D[0] = 1e-10
        sigma2 = (D2-Dm**2)
        mu = Dm
        pi = math.acos(-1.0)
        dsd = np.zeros(nd)
        if (Dm > 0.0 and D2 > 0.0):
                B = 1.0/np.sqrt(2.0*pi*sigma2)
                f = B*np.exp(-(D-mu)**2/2.0/sigma2)
                norm = sum(f*h)
                if (norm > 0.0) :
                        dsd[:] = f/norm
        return [D, dsd]


Dm = 1.0
D2 = 0.01

[D,dsd] = computeDSDrr(Dm,D2);


plt.plot(D,dsd)
plt.show()
