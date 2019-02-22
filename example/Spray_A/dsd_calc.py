import numpy as np
import matplotlib.pyplot as plt
import csv

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
with open('dsd.csv', 'w') as f:
    i = 0
    for val in di:
        f.write("%1.16E,%1.16E\n" % val,md[i])
	i = i + 1
'''
hist, bin_edges = np.histogram(np.array(di)/Dn,bins=100000,density=True)
 
#plt.hist(np.array(md), bins=[0,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1])
#plt.xlim([0,4])
plt.plot(bin_edges[:-1],hist)
plt.show()
