import sys
from scipy.special import *
import math
import numpy as np
from matplotlib import pyplot as plt

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

fsz=24
plt.rc('text', usetex=True)
font = {'family':'serif','size':fsz, 'serif': ['Computer Modern Roman']}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':fsz})

def computeDSDgg3(Dm,D2,D3):
        L = 10.0*Dm
        nd = 200
        h = L/nd
        D = h*np.arange(0,nd)
        D[0] = 1e-10
	var = D2-Dm**2
        stdd = np.sqrt(var)
        skew = (D3-3*Dm*D2+2*Dm**3)/var**1.5
        #print "Mean: ",Dm, "Var: ",var, "Skewness: ",skew
        #[d,p] = bisection2(ratio1,ratio2)
        mu = Dm
	sigma = stdd
        Q = skew
	Qm2 = Q**(-2)
	w = np.log(Q**2*gamma(Qm2))/Q
	B = abs(Q)*(Qm2)**(Qm2)/sigma/gamma(Qm2)*np.exp(Qm2*(Q*w-np.exp(Q*w)))
        dsd = np.zeros(nd)
	if (Dm > 0.0 and D2 > 0.0):
	        f = B/D
        	norm = sum(f*h)
        	if (norm > 0.0):
                	dsd[:] = f/norm
        return [D, dsd]

def computeDSDgg(Dm,D2,D3):
        L = 10.0*Dm
        nd = 200
        h = L/nd
        D = h*np.arange(0,nd)
        D[0] = 1e-10
        var = D2-Dm**2
	stdd = np.sqrt(var)
	skew = (D3-3*Dm*D2+2*Dm**3)/var**1.5
	#print "Mean: ",Dm, "Var: ",var, "Skewness: ",skew
        #[d,p] = bisection2(ratio1,ratio2)
	shape = 4.0/skew**2
	scale = stdd/np.sqrt(shape)
	shift = Dm - 2.0*stdd/skew
	alpha = shape
	beta = scale
	gamm = shift
	B = 1.0/(beta**alpha)/gamma(alpha)
	print alpha, beta, gamm, B
        dsd = np.zeros(nd)
        if (Dm > 0.0 and D2 > 0.0):
                f = B*((D-gamm)**(alpha-1.0))*(np.exp(-((D-gamm)/beta)))
		norm = sum(f*h)
                if (norm > 0.0):
                        dsd[:] = f/norm
	dsd[0] = 0.0
	return [D, dsd]

def computeDSDpn(Dm,D2,D3):
        L = 10*Dm
        nd = 200
        h = L/nd
        D = h*np.arange(0,nd)
        D[0] = 1e-10
        lamda = Dm
	var = D2-Dm**2
	skew = (D3-3*Dm*D2+2*Dm**3)/var**1.5
	#print "Var: ",var, "Skewness: ",skew
        alpha = np.sqrt(skew*var)
	m = 0.5*(skew+1.0)
        B = 1.0/alpha/beta(skew/2,0.5)
        dsd = np.zeros(nd)
        if (Dm > 0.0 and D2 > 0.0):
                f = B*(1.0+((D-Dm)/alpha)**2)**(-0.5*(skew+1.0))
                norm = sum(f*h)
                if (norm > 0.0) :
                        dsd[:] = f/norm
        return [D, dsd]

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

def computeDSD(Dm,D2):
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

def computeDSDln(Dm,D2):
        L = 10*Dm
        nd = 200
        h = L/nd
        D = h*np.arange(0,nd)
	D[0] = 1e-10
	mu = np.log(Dm**2/np.sqrt(D2))
	sigma2 = np.log(D2/Dm**2)
	pi = math.acos(-1.0)
	B = 1.0/np.sqrt(2.0*pi*sigma2)
        dsd = np.zeros(nd)
        if (Dm > 0.0 and D2 > 0.0):
                f = B*np.exp(-((np.log(D)-mu)**2/2/sigma2))/D 
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

def computeDSDgd(Dm,D2):
        L = 10*Dm
        nd = 200
        h = L/nd
        D = h*np.arange(0,nd)
        D[0] = 1e-10
        alpha = Dm**2/(D2-Dm**2)
        beta = Dm/(D2-Dm**2)
        pi = math.acos(-1.0)
        dsd = np.zeros(nd)
        if (Dm > 0.0 and D2 > 0.0):
        	B = beta**alpha/gamma(alpha)
	        f = B*D**(alpha-1.0)*np.exp(-beta*D)
                norm = sum(f*h)
                if (norm > 0.0) :
                        dsd[:] = f/norm
        return [D, dsd]

nozD = 90e-6
dataRR = np.loadtxt(sys.argv[1],comments="x", unpack=False)
dd21 = np.loadtxt("DSD_1mm.dat",comments="#", unpack=False)
dd22 = np.loadtxt("DSD_2mm.dat",comments="#", unpack=False)
dd24 = np.loadtxt("DSD_4mm.dat",comments="#", unpack=False)
dd26 = np.loadtxt("DSD_6mm.dat",comments="#", unpack=False)
dd28 = np.loadtxt("DSD_8mm.dat",comments="#", unpack=False)
dd210 = np.loadtxt("DSD_10mm.dat",comments="#", unpack=False)
dd211 = np.loadtxt("DSD_11mm.dat",comments="#", unpack=False)

dd21[:,0] = dd21[:,0]*nozD

nez=np.max(np.shape(dd21))
suum = 0.0
for i in np.arange(0,nez-1):
        hh = dd21[:,0]/nozD
        ff = dd21[:,1]
        suum = suum + 0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

dd21[:,0] = hh
dd21[:,1] = ff/suum

dm21 = 0.0
for i in np.arange(0,nez-1):
        hh = dd21[:,0]
        ff = dd21[:,1]
        dm21 = dm21 + 0.5*(hh[i]+hh[i+1])*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d221 = 0.0
for i in np.arange(0,nez-1):
        hh = dd21[:,0]
        ff = dd21[:,1]
        d221 = d221 + (0.25*(hh[i]+hh[i+1])**2)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d321 = 0.0
for i in np.arange(0,nez-1):
        hh = dd21[:,0]
        ff = dd21[:,1]
        d321 = d321 + ((0.5*(hh[i]+hh[i+1]))**3)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

print dm21, d221, d321

dd22[:,0] = dd22[:,0]

nez=np.max(np.shape(dd22))
suum = 0.0
for i in np.arange(0,nez-1):
	hh = dd22[:,0]/nozD
	ff = dd22[:,1]
	suum = suum + 0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])
	
dd22[:,0] = hh
dd22[:,1] = ff/suum

dm22 = 0.0
for i in np.arange(0,nez-1):
        hh = dd22[:,0]
        ff = dd22[:,1]
        dm22 = dm22 + 0.5*(hh[i]+hh[i+1])*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d222 = 0.0
for i in np.arange(0,nez-1):
        hh = dd22[:,0]
        ff = dd22[:,1]
        d222 = d222 + (0.25*(hh[i]+hh[i+1])**2)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d322 = 0.0
for i in np.arange(0,nez-1):
        hh = dd22[:,0]
        ff = dd22[:,1]
        d322 = d322 + ((0.5*(hh[i]+hh[i+1]))**3)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

print dm22, d222, d322

[D22,DSD22] = computeDSD(dm22,d222)
[D22l,DSD22l] = computeDSDln(dm22,d222)
[D22rl,DSD22rl] = computeDSDrl(dm22,d222)
[D22gd,DSD22gd] = computeDSDgd(dm22,d222)
[D22gg,DSD22gg] = computeDSDgg(dm22,d222,d322)

nez=np.max(np.shape(dd24))
suum = 0.0
for i in np.arange(0,nez-1):
        hh = dd24[:,0]/nozD
        ff = dd24[:,1]
        suum = suum + 0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

dd24[:,0] = hh
dd24[:,1] = ff/suum

dm24 = 0.0
for i in np.arange(0,nez-1):
        hh = dd24[:,0]
        ff = dd24[:,1]
        dm24 = dm24 + 0.5*(hh[i]+hh[i+1])*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d224 = 0.0
for i in np.arange(0,nez-1):
        hh = dd24[:,0]
        ff = dd24[:,1]
        d224 = d224 + (0.25*(hh[i]+hh[i+1])**2)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d324 = 0.0
for i in np.arange(0,nez-1):
        hh = dd24[:,0]
        ff = dd24[:,1]
        d324 = d324 + ((0.5*(hh[i]+hh[i+1]))**3)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

[D24,DSD24] = computeDSD(dm24,d224)
[D24l,DSD24l] = computeDSDln(dm24,d224)
[D24rl,DSD24rl] = computeDSDrl(dm24,d224)
[D24gd,DSD24gd] = computeDSDgd(dm24,d224)
[D24gg,DSD24gg] = computeDSDgg(dm24,d224,d324)

nez=np.max(np.shape(dd26))
suum = 0.0
for i in np.arange(0,nez-1):
        hh = dd26[:,0]/nozD
        ff = dd26[:,1]
        suum = suum + 0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

dd26[:,0] = hh
dd26[:,1] = ff/suum

dm26 = 0.0
for i in np.arange(0,nez-1):
        hh = dd26[:,0]
        ff = dd26[:,1]
        dm26 = dm26 + 0.5*(hh[i]+hh[i+1])*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d226 = 0.0
for i in np.arange(0,nez-1):
        hh = dd26[:,0]
        ff = dd26[:,1]
        d226 = d226 + (0.25*(hh[i]+hh[i+1])**2)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d326 = 0.0
for i in np.arange(0,nez-1):
        hh = dd26[:,0]
        ff = dd26[:,1]
        d326 = d326 + ((0.5*(hh[i]+hh[i+1]))**3)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

[D26,DSD26] = computeDSD(dm26,d226)
[D26l,DSD26l] = computeDSDln(dm26,d226)
[D26rl,DSD26rl] = computeDSDrl(dm26,d226)
[D26gd,DSD26gd] = computeDSDgd(dm26,d226)
[D26gg,DSD26gg] = computeDSDgg(dm26,d226,d326)

nez=np.max(np.shape(dd28))
suum = 0.0
for i in np.arange(0,nez-1):
        hh = dd28[:,0]/nozD
        ff = dd28[:,1]
        suum = suum + 0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

dd28[:,0] = hh
dd28[:,1] = ff/suum

dm28 = 0.0
for i in np.arange(0,nez-1):
        hh = dd28[:,0]
        ff = dd28[:,1]
        dm28 = dm28 + 0.5*(hh[i]+hh[i+1])*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d228 = 0.0
for i in np.arange(0,nez-1):
        hh = dd28[:,0]
        ff = dd28[:,1]
        d228 = d228 + (0.25*(hh[i]+hh[i+1])**2)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d328 = 0.0
for i in np.arange(0,nez-1):
        hh = dd28[:,0]
        ff = dd28[:,1]
        d328 = d328 + ((0.5*(hh[i]+hh[i+1]))**3)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

[D28,DSD28] = computeDSD(dm28,d228)
[D28l,DSD28l] = computeDSDln(dm28,d228)
[D28rl,DSD28rl] = computeDSDrl(dm28,d228)
[D28gd,DSD28gd] = computeDSDgd(dm28,d228)
[D28gg,DSD28gg] = computeDSDgg(dm28,d228,d328)

nez=np.max(np.shape(dd210))
suum = 0.0
for i in np.arange(0,nez-1):
        hh = dd210[:,0]/nozD
        ff = dd210[:,1]
        suum = suum + 0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

dd210[:,0] = hh
dd210[:,1] = ff/suum

dm210 = 0.0
for i in np.arange(0,nez-1):
        hh = dd210[:,0]
        ff = dd210[:,1]
        dm210 = dm210 + 0.5*(hh[i]+hh[i+1])*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d2210 = 0.0
for i in np.arange(0,nez-1):
        hh = dd210[:,0]
        ff = dd210[:,1]
        d2210 = d2210 + (0.25*(hh[i]+hh[i+1])**2)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d3210 = 0.0
for i in np.arange(0,nez-1):
        hh = dd210[:,0]
        ff = dd210[:,1]
        d3210 = d3210 + ((0.5*(hh[i]+hh[i+1]))**3)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

[D210,DSD210] = computeDSD(dm210,d2210)
[D210l,DSD210l] = computeDSDln(dm210,d2210)
[D210rl,DSD210rl] = computeDSDrl(dm210,d2210)
[D210gd,DSD210gd] = computeDSDgd(dm210,d2210)
[D210gg,DSD210gg] = computeDSDgg(dm210,d2210,d3210)

nez=np.max(np.shape(dd211))
suum = 0.0
for i in np.arange(0,nez-1):
        hh = dd211[:,0]/nozD
        ff = dd211[:,1]
        suum = suum + 0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

dd211[:,0] = hh
dd211[:,1] = ff/suum

dm211 = 0.0
for i in np.arange(0,nez-1):
        hh = dd211[:,0]
        ff = dd211[:,1]
        dm211 = dm211 + 0.5*(hh[i]+hh[i+1])*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d2211 = 0.0
for i in np.arange(0,nez-1):
        hh = dd211[:,0]
        ff = dd211[:,1]
        d2211 = d2211 + (0.25*(hh[i]+hh[i+1])**2)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

d3211 = 0.0
for i in np.arange(0,nez-1):
        hh = dd211[:,0]
        ff = dd211[:,1]
        d3211 = d3211 + ((0.5*(hh[i]+hh[i+1]))**3)*0.5*(ff[i]+ff[i+1])*(hh[i+1]-hh[i])

[D211,DSD211] = computeDSD(dm211,d2211)
[D211l,DSD211l] = computeDSDln(dm211,d2211)
[D211rl,DSD211rl] = computeDSDrl(dm211,d2211)
[D211gd,DSD211gd] = computeDSDgd(dm211,d2211)
[D211gg,DSD211gg] = computeDSDgg(dm211,d2211,d3211)

SOI = 0.0

print dataRR[np.round(2.0e-3/nozD),0],dm22, d222
print dataRR[np.round(4.0e-3/nozD),0],dm24, d224
print dataRR[np.round(6.0e-3/nozD),0],dm26, d226
print dataRR[np.round(8.0e-3/nozD),0],dm28, d228
print dataRR[np.round(10.0e-3/nozD),0],dm210, d2210
print dataRR[np.round(11.0e-3/nozD),0],dm211, d2211

# Droplet Size Distribution
nd = 200
nx = max(np.shape(dataRR[:,0]))
RR22 = np.zeros([nd,2])
RR24 = np.zeros([nd,2])
RR26 = np.zeros([nd,2])
RR28 = np.zeros([nd,2])
RR210 = np.zeros([nd,2])
Z = dataRR[:,0]
Dm = dataRR[:,8]
D2 = dataRR[:,9]
for k in np.arange(0,nx):
        if (Dm[k] > 0.0 and D2[k] > 0.0):
                if (Z[k] == np.round(2.0e-3/nozD)):
                        [RR22[:,0],RR22[:,1]] = computeDSD(Dm[k],D2[k])
                if (Z[k] == np.round(4.0e-3/nozD)):
                        [RR24[:,0],RR24[:,1]] = computeDSD(Dm[k],D2[k])
                if (Z[k] == np.round(6.0e-3/nozD)):
                        [RR26[:,0],RR26[:,1]] = computeDSD(Dm[k],D2[k])
                if (Z[k] == np.round(8.0e-3/nozD)):
                        [RR28[:,0],RR28[:,1]] = computeDSD(Dm[k],D2[k])
                if (Z[k] == np.round(10.0e-3/nozD)):
                        [RR210[:,0],RR210[:,1]] = computeDSD(Dm[k],D2[k])

plt.plot(RR22[:,0]*nozD*1e6,RR22[:,1]/(nozD*1e6), color='red', ls='-', lw=2, label='Sim: Rosin-Rammler')

plt.plot(dd22[:,0]*nozD*1e6,dd22[:,1]/(nozD*1e6), linestyle="None",marker='s',mfc='None',mec='red',mew=2,markevery=1, label='Reference')
plt.ylabel('Distribution [1/$\mu m$]')
plt.xlabel('d [$\mu m$]')
plt.xlim([0,4])
#plt.ylim([0,120])
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('DSD_2mm.pdf', bbox_inches='tight')
plt.show()

plt.clf()
plt.plot(RR24[:,0]*nozD*1e6,RR24[:,1]/(nozD*1e6), color='green', ls='-', lw=2, label='Sim: Rosin-Rammler')
plt.plot(dd24[:,0]*nozD*1e6,dd24[:,1]/(nozD*1e6), linestyle="None",marker='s',mfc='None',mec='green',mew=2,markevery=1, label='Reference')
plt.ylabel('Distribution [1/$\mu m$]')
plt.xlabel('d [$\mu m$]')
plt.xlim([0,4])
#plt.ylim([0,120])
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('DSD_4mm.pdf', bbox_inches='tight')
plt.show()

plt.clf()
plt.plot(RR26[:,0]*nozD*1e6,RR26[:,1]/(nozD*1e6), color='blue', ls='-', lw=2, label='Sim: Rosin-Rammler')
plt.plot(dd26[:,0]*nozD*1e6,dd26[:,1]/(nozD*1e6), linestyle="None",marker='s',mfc='None',mec='blue',mew=2,markevery=1, label='Reference')
plt.ylabel('Distribution [1/$\mu m$]')
plt.xlabel('d [$\mu m$]')
plt.xlim([0,4])
#plt.ylim([0,120])
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('DSD_6mm.pdf', bbox_inches='tight')
plt.show()

plt.clf()
plt.plot(RR28[:,0]*nozD*1e6,RR28[:,1]/(nozD*1e6), color='magenta', ls='-', lw=2, label='Sim: Rosin-Rammler')
plt.plot(dd28[:,0]*nozD*1e6,dd28[:,1]/(nozD*1e6), linestyle="None",marker='s',mfc='None',mec='magenta',mew=2,markevery=1, label='Reference')
plt.ylabel('Distribution [1/$\mu m$]')
plt.xlabel('d [$\mu m$]')
plt.xlim([0,4])
#plt.ylim([0,120])
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('DSD_8mm.pdf', bbox_inches='tight')
plt.show()

plt.clf()
plt.plot(RR210[:,0]*nozD*1e6,RR210[:,1]/(nozD*1e6), color='purple', ls='-', lw=2, label='Sim: Rosin-Rammler')
plt.plot(dd210[:,0]*nozD*1e6,dd210[:,1]/(nozD*1e6), linestyle="None",marker='s',mfc='None',mec='purple',mew=2,markevery=1, label='Reference')
plt.ylabel('Distribution [1/$\mu m$]')
plt.xlabel('d [$\mu m$]')
plt.xlim([0,4])
#plt.ylim([0,120])
legend = plt.legend(loc='upper right', shadow=True)
plt.savefig('DSD_10mm.pdf', bbox_inches='tight')
plt.show()
