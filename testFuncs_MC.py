import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt

#import IR_LightEchoes_NewMeth as IRLE
from FluxFuncs_IRLE import *



Plot = False
ResStudy = False

nne = 1.
nu0 = numicron

#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
Omb = 1./(1*yr2sec)
L0 = 6.78*10**46
MPGmx = 10**9.4*Msun
Ryr = c*yr2sec
RdPG = np.sqrt(0.1)*2.8 *pc2cm
OmPG = Omb*2.*np.pi/4.1
alphnu = 1.1

Rorb = c*2.*np.pi/Omb
Ompc = 2.*np.pi*c/pc2cm/2.


## TEST VALUES
Lav = L0
betst = 0.10
Inc = ma.acos(0.067/betst)#0.*np.pi/4.
Ombn = OmPG
alph = 0.0

Rde = RdPG
pp = 2.0
thetTst = 1.*np.pi/4
JJt =np.pi/4
aeff = 0.16*10**(-4) #(0.1 micrometer is an average ISM dust grain size)


Dst = 1.4*10**9*pc2cm
Rrout = 2.0*Rde

md = 10**(-14)
n0 = 1.2/(ma.pi * aeff*aeff * Rde)



W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3

nuVbnd = c/(545*10**(-7))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-20)*(W1mn + W1mx)/2
FW2Rel = 1.7187*10**(-20)*(W2mn + W2mx)/2


Nt=20
Nr=10
Nth=100
Nph=100
tt = np.linspace(0., 2.,       Nt)*2*np.pi/Ombn
rr = np.linspace(1., 10.,      Nr)*Rde
th = np.linspace(0., np.pi,    Nth)
ph = np.linspace(0., 2.*np.pi, Nph)
nnu = np.linspace(0.01, 2.,       Nt)*numicron

##TABULATE T's and RHSs
print "Creating look up tables"
NT = 1000
NMC = 1.e5
RHS_table = np.zeros(NT)
T_table = np.linspace(1., 2000., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)



Dop_args = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]

Amp = 0.14
t0 = 0.0
ISO_args = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]


###-------------profiling-----------####
# Nt=10
# for i in range (0, Nt):
# 	#Fnuplt[i] = -2.5*np.log10(Fnudphidth(Rde, numicron, tt[i], Dst, Rrout, Targs)/FW1Rel)
# 	Fnu(numicron, 0., Dst, Rrout, Targs)

# import sys
# sys.exit(0)
###-------------profiling-----------####
print "Plotting stuff"


Nx = 100
xx = np.linspace(-1.*Rde, 1.*Rde, Nx)
zz = np.linspace(-1.*Rde, 1.*Rde, Nx)

#xx = Rde* np.sin(th) * np.cos(ph)
#zz = Rde* np.cos(th)
#ytst = ( Rde*Rde - (xx*xx + zz*zz) )**(0.5)

ytst = 0*Rde

if (Plot):
	tauT = np.zeros([Nx,Nx])
	nd = np.zeros([Nx,Nx])
	for i in range(0, Nx):
	 	for j in range(0, Nx):
			#ytst = ( Rde*Rde - (xx[i]*xx[i] + zz[j]*zz[j]) )**(0.5)
			#tauObs(nu, x, y, z, Rout, aeff, n0, Rd, p, thetT, JJ, nu0, nn)
	 		#tauT[j][i] = tauObs_Shell(2.*numicron, xx[i], ytst, zz[j], Rrout, aeff, n0, Rde, 2., thetTst, JJt, numicron, 1.)
	 		nd[j][i] =  np.log10(nDust(xx[i], ytst, zz[j], n0, Rde, 2.0, thetTst, JJt))
	#nd =  nDust(xx, ytst, zz, n0, Rde, 2.0, thetTst, JJt)

	TD_angles = np.empty([Nth,Nph])
	TD_rth = np.empty([Nth,Nr])
	for i in range(0, Nth):
		for j in range(0, Nph):
			#TDust_Dop(t,r,thet,phi,args, RHStable, Ttable)
			TD_angles[i][j] =  TDust_Dop(0.5*2.*np.pi/Ombn, Rde, th[i], ph[j], Dop_args, RHS_table, T_table) 



	for i in range(0, Nth):
		for j in range(0, Nr):
			TD_rth[i][j] =  TDust_Dop(0.5*2.*np.pi/Ombn, rr[j], th[i], 0, Dop_args, RHS_table, T_table) 



 	FobsThick = np.zeros(Nt)
	for i in range(Nt):
	 	FobsThick[i] = -2.5*np.log10(F_Thick_Dop_QuadInt(W1mn, W1mx, tt[i], Dst, Dop_args, RHS_table, T_table)/FW1Rel)
	# 	FobsThick[i] = -2.5*np.log10(Thick_Dop_MCInt(NMC, W1mn, W1mx, tt[i], Dst, Dop_args, RHS_table, T_table)/FW1Rel)



	plt.figure()
	plt.plot(tt/(2*np.pi/Ombn),FobsThick)
	plt.show()

	plt.figure()
	#plt.subplot(211)
	#tauplt = plt.contourf(xx, zz,  tauT, colorbar=True ) 
	#plt.colorbar(tauplt)
	#plt.subplot(212)
	ndplt = plt.contourf(xx, zz,  nd, colorbar=True ) 
	plt.colorbar(ndplt)
	plt.show()


	plt.figure()
	plt.subplot(121)
	Ta = plt.contourf(ph,th,TD_angles, contours=100, colorbar = True)
	plt.colorbar(Ta)
	plt.subplot(122)
	Tr = plt.contourf(rr,th,TD_rth, contours=100, colorbar = True)
	plt.colorbar(Tr)
	plt.show()


import time
print "start timing"


t1=time.clock()
#chk1=-2.5*np.log10(Thick_Dop_MCInt(NMC, W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)/FW1Rel)
#chk1 = Thick_Dop_MCInt(NMC, W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)
t2=time.clock()
td = (t2-t1)/2.
print "Fobs_Thick MCint time = %g" %td
#print "Fobs_Thick MCint Value = %g" %chk1

t1=time.clock()
chk2=-2.5*np.log10(F_Thick_Dop_QuadInt(W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)/FW1Rel)
#chk2 = F_Thick_Dop_QuadInt(W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)
t2=time.clock()
td = (t2-t1)/2.
print "Fobs_Thick Quad time = %g" %td
#print "Fobs_Thick Quad Value = %g" %chk2

if (ResStudy):
	Nres = 10
	res = []
	NMC = np.logspace(5,6,Nres)
	for i in range(Nres):
		chk=-2.5*np.log10(Thick_Dop_MCInt(NMC[i], W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)/FW1Rel)
		#chk=Thick_Dop_MCInt(NMC[i], W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)
		res.append(chk)


	plt.figure()
	plt.plot(NMC, res)
	plt.show()

# Nres = 10
# res = []
# for i in range(Nres):
# 	#chk2=-2.5*np.log10(F_Thick_Dop_QuadInt(W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)/FW1Rel)
# 	chk0=-2.5*np.log10(Thick_Dop_MCInt(1.e3, W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)/FW1Rel)
# 	chk1=-2.5*np.log10(Thick_Dop_MCInt(1.e4, W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)/FW1Rel)
# 	chk2=-2.5*np.log10(Thick_Dop_MCInt(1.e5, W1mn, W1mx, 0.0, Dst, Dop_args, RHS_table, T_table)/FW1Rel)
# 	res1.append(chk1)
# 	res2.append(chk2)
# 	res3.append(chk3)

# mn1 = np.mean(res1)
# sd1 = np.std(res1)

# mn2 = np.mean(res2)
# sd2 = np.std(res2)

# mn3 = np.mean(res3)
# sd3 = np.std(res3)









