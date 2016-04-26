import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt

import IR_LightEchoes_NewMeth as IRLE
from IR_LightEchoes_NewMeth import *


nne = 1.
nu0 = numicron*0.2

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
Inc = ma.acos(0.07/betst)#0.*np.pi/4.
Ombn = OmPG
alph = 1.1

Rde = RdPG
pp = 2.0
thetTst = 0.*np.pi/4
JJt =-np.pi/8
aeff = 0.1*10**(-4) #(0.1 micrometer is an average ISM dust grain size)


Dst = 1.4*10**9*pc2cm
Rrout = 1.0*Rde

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
RHS_table = np.zeros(NT)
T_table = np.linspace(1., 2000., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)

print "compute tau_los tables"
nx=101
nz=101
xx = np.linspace(-Rrout, Rrout, nx)
zz = np.linspace(-Rrout, Rrout, nz)
#ph = np.linspace(0, 2.*ma.pi, 10)
#rrg,thg,phg = np.meshgrid(rr, th, ph)


# tauGrid = [[],[],[]]
# for i in range(0, nx):
# 	for j in range(0, nz):
# 			tauGrid[0].append(tauObs(xx[i], zz[j], Rrout, aeff, n0, Rde, pp, thetTst, JJt))
# 			tauGrid[1].append(xx[i])
# 			tauGrid[2].append(zz[j])




TD_angles = np.empty([Nth,Nph])
TD_rth = np.empty([Nth,Nr])
FsrcPlt = np.empty(Nt)
FnuW1 = np.empty(Nt)
FnuW2 = np.empty(Nt)
Fnuplt = np.empty(Nt)

Targs = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]


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


tauT = np.zeros([Nx,Nx])
nd = np.zeros([Nx,Nx])
for i in range(0, Nx):
	for j in range(0, Nx):
		#ytst = ( Rde*Rde - (xx[i]*xx[i] + zz[j]*zz[j]) )**(0.5)
		#tauObs(nu, x, y, z, Rout, aeff, n0, Rd, p, thetT, JJ, nu0, nn)
 		tauT[j][i] = tauObs_Shell(2.*numicron, xx[i], ytst, zz[j], Rrout, aeff, n0, Rde, 2., thetTst, JJt, numicron, 1.)
 		nd[j][i] =  np.log10(nDust(xx[i], ytst, zz[j], n0, Rde, 2.0, thetTst, JJt))




plt.figure()
plt.subplot(211)
tauplt = plt.contourf(xx, zz,  tauT, colorbar=True ) 
plt.colorbar(tauplt)
plt.subplot(212)
ndplt = plt.contourf(xx, zz,  nd, colorbar=True ) 
plt.colorbar(ndplt)
plt.show()

# for i in range(0, len(tt)):
# 	Targs_t.append([tt[i], Rde, np.pi/2., 0.0, Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne])
# 	TD_t.append( IRLE.TDust(Targs_t[i]) )


for i in range(0, Nth):
	for j in range(0, Nph):
		#Targs_angles = [0., Rde, th[i], ph[j], Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		TD_angles[i][j] =  IRLE.TDust(0.5*2.*np.pi/Ombn, 2.1*Rde, th[i], ph[j], Targs, RHS_table, T_table) 



for i in range(0, Nth):
	for j in range(0, Nr):
#		Targs_rth = [0., rr[j], th[i], 0., Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		TD_rth[i][j] =  IRLE.TDust(0.5*2.*np.pi/Ombn, rr[j], th[i], 0, Targs, RHS_table, T_table) 



# for i in range (0, Nt):
# # 	#Fnuplt[i] = -2.5*np.log10(Fnu(numicron, tt[i], Dst, Rrout, Targs)/FW1Rel)
# # 	#Fnuplt[i] = -2.5*np.log10(Fnudphidth(Rde, numicron, tt[i], Dst, Rrout, Targs)/FW1Rel)
# 	FsrcPlt[i] = -2.5*np.log10(Fsrc(tt[i], Dst, ma.pi/2., 0.0, Lav/14., betst, Inc, Ombn, alph)/FVbndRel)-3.0
#  	FnuW1[i] = -2.5*np.log10(Fobs_Shell(W1mn, W1mx, tt[i], Dst, Rrout, Targs, RHS_table, T_table)/FW1Rel)-0.0
#  	FnuW2[i] = -2.5*np.log10(Fobs_Shell(W2mn, W2mx, tt[i], Dst, Rrout, Targs, RHS_table, T_table)/FW2Rel)+0.5
# # 	#FnuW1[i] = -2.5*np.log10(FobsR(W1mn, W1mx, tt[i], Dst, Rrout, Targs)/FW1Rel)


# plt.figure()
# plt.plot(tt/(2*np.pi/Omb), FsrcPlt)
# plt.plot(tt/(2*np.pi/Omb), FnuW1)
# plt.plot(tt/(2*np.pi/Omb), FnuW2)
# plt.show()


#TD_angles =  IRLE.TDust(0., Rde, th, ph, Targs) 


#Targs_angles = [ Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
#TD_angles =  IRLE.TDust(0., rr, th, 0,Targs_rth) 

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

xt = 4*Rde
zt = 0.0
thet = ma.pi/2
ph = 0.0
r = ma.sqrt(xt*xt + zt*zt)
t1=time.clock()

xe     = Rrout*( 1. - (r/Rrout)*(r/Rrout) * (  np.cos(thet)*np.cos(thet)  +  np.sin(thet)*np.sin(ph) * np.sin(thet)*np.sin(ph)  )  )**(0.5)
tauObs = np.pi*aeff*aeff * intg.quad(nDust  ,xt, xe , args=(0., zt, n0, Rde, pp, thetTst, JJt) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

# ixl = np.where( xt < tauGrid[1])[0].max()
# ixall = np.where(tauGrid[1][ixl] ==  tauGrid[1])[0]

# TG = np.transpose(tauGrid)
# TGz = np.transpose(TG[ixall])  ## this is all of the z values for given x value
# it =np.where(zt < TGz[2])[0].max()

# tauObs = tauGrid[0][it]	
t2=time.clock()
print t2-t1


t1=time.clock()
chk0=Fnuint_Shell(0.0, ma.pi/2., W1mn, Rde,1.0, Dst, Rrout, Targs, RHS_table, T_table)
t2=time.clock()
td = t2-t1
print "Fnuint_Shell time = %g" %td

t1=time.clock()
chk1=Fnudphi_Shell(ma.pi/2, W1mn, Rde,0.0, Dst, Rrout, Targs, RHS_table, T_table)
t2=time.clock()
td = t2-t1
print "Fnudphi_Shell time = %g" %td

t1=time.clock()
chk2=Fnu_Shell(W1mn, Rde,0.0, Dst, Rrout, Targs, RHS_table, T_table)
t2=time.clock()
td = t2-t1
print "Fnu_Shell time = %g" %td


t1=time.clock()
chk3=-2.5*np.log10(Fobs_Shell(W1mn,W1mx, Rde,[0.,1.], Dst, Rrout, Targs, RHS_table, T_table)/FW1Rel)
t2=time.clock()
td = (t2-t1)/2.
print "Fobs_Shell timw= %g" %td


t1=time.clock()
chk4=-2.5*np.log10(Fobs_Thick(W1mn,W1mx, [1.,2.], Dst, Rrout, Targs, RHS_table, T_table)/FW1Rel)
t2=time.clock()
td = (t2-t1)/2.
print "Fobs_Thick timw= %g" %td




