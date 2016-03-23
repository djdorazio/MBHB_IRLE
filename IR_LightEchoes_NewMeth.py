import numpy as np
import math as ma
import numexpr as ne
import scipy as sc
#from scipy.special import gamma as Gamma
#from scipy.special import expn as ExpIntegralE  ## exponential integral (z,n) same as mathematica
from scipy.optimize import brentq #fmin

import scipy.integrate as intg
#from apw import Tfunc

#Goal is to out put the IR flux from irradiated dust geometry surrouding binary.
#Result Function 1: gives nu and r dependnet flux
#Result Function 2: gives flux integrated over r and nu - given bounds

#### INTEGRATION ERROR TOLS
myrel = 0.1
myabs = 0.1
reclim = 1
limlst = 3
maxp1 = 1
fo = 1

##GLOBAL PHYSICS CONSTANTS (cgs):
c = 2.9979*10**(10)
G = 6.673*10**(-8)
Msun = 1.998*10**(33)

kb = 1.3807*10**(-16)
h = 6.62607*10**(-27) 
sigSB = 5.670*10**(-5)

pc2cm = 3.08567758*10**(18)
numicron = c/(10**(-4))

yr2sec = 3600.*24.*365.25



## Back body specific intensity
def Bv(nu, T):
	return 2.*h*nu*nu*nu/(c*c)*1./(np.exp( h*nu/(kb*T) ) - 1.)
## Dust absorption efficiency
def Qv(nu, nu0, nn):
	#qv = (nu/nu0)**(nn)
	#qv = min(qv, 1.0)
	#if (type(qv) is float):
	#	if (qv>1.0):
	#		qv=1.0
	#else:
	#	ii = np.where(qv > 1.0)[0]
	#	qv[ii] = 1.0
	return np.min( [(nu/nu0)**(nn), 1.*nu/nu])#qv #
	


def QvBv(nu, T, nu0, nn):
	qv = (nu/nu0)**(nn)
	if (type(qv) is float):
		qv = min(qv, 1.0)
		#if (qv>1.0):
		#	qv=1.0
	else:
		ii = np.where(qv > 1.0)[0]
		qv[ii] = 1.0
	return 2.*h*nu*nu*nu/(c*c)*1./(np.exp(h*nu/(kb*T)) - 1.) * qv



# Torical dust profile
def nDust(x,y,z, n0, Rd, p, thetT, JJ):
	nprof = 0.0
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	rofx  = (x*x + y*y + z*z)**(0.5)
	throt = np.arctan2(np.sqrt(xrot*xrot + y*y), zrot)
	if (rofx>=Rd and throt>thetT and throt<(np.pi - thetT)):
		nprof = n0*(rofx/Rd)**(-p)

	return nprof


## flux from beaming binary as seen by observer at r, theta, phi (centere on binary barycenter)
def Fsrc(t, r, thet, phi, Lavg, bets, incl, Ombin, alphnu):
#Rot by Inc around y axis, Rotation around Bin ang momentum axis by Ombin*t 
#starting point of secindary at t=0
	phis = -np.pi/2
	thetas = np.pi/2

	Vxorb = np.cos(incl) *np.sin(phis + Ombin*t) * np.sin(thetas)
	Vyorb = -np.cos(phis + Ombin*t) * np.sin(thetas)
	Vzorb = -np.sin(incl)* np.sin(phis + Ombin*t) * np.sin(thetas)

#Unit Position of observer on dust shell in dust shell coords
	xobs = np.sin(thet)*np.cos(phi)
	yobs = np.sin(thet)*np.sin(phi)
	zobs = np.cos(thet)

#Line of sight velocity anywhere in dust shell (vorb.rhatobs)	
	FracLOS = Vxorb*xobs + Vyorb*yobs + Vzorb*zobs

#Doppler factor 
	Gams = (1. - bets*bets)**(-0.5)
	Dop = 1./(Gams*(1. - bets*FracLOS))

# return flux at observer coordinates r, phi, theta (in dust or at earth)
	return Lavg/(4.*np.pi*r*r)*(Dop)**(3. - alphnu)


## equation to tabulate RHS and T
def T_RHS(Td, nu0, nn):
	RHS = 4.*np.pi*  (intg.quad(QvBv  ,0., nu0 , args=(Td, nu0, nn) )[0] + intg.quad(Bv  ,nu0 ,numicron*1000., args=(Td) )[0])#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	#return np.log10(RHS)
	return RHS



	


def TDust(t,r,thet,phi,args, RHStable, Ttable):
	#thetT = args[8] 
	#JJ = args[9] 
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	x = r*np.sin(thet)*np.cos(phi)
	y = r*np.sin(thet)*np.sin(phi)
	z = r*np.cos(thet)
	rofx  = (x*x + y*y + z*z)**(0.5) 
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	Tprof = 0.01*t/t
	if (rofx>=Rd and throt>thetT and throt<(np.pi - thetT)):

	###-----------------###
	### COMPUTE Fsrc    ###
	###-----------------###
		#Rot by Inc around y axis, Rotation around Bin ang momentum axis by Ombin*t 
		#starting point of secindary at t=0
		phis = -np.pi/2.
		thetas = np.pi/2.


		Vxorb = np.cos(incl) *np.sin(phis + Ombin*t) * np.sin(thetas)
		Vyorb = -np.cos(phis + Ombin*t) * np.sin(thetas)
		Vzorb = -np.sin(incl)* np.sin(phis + Ombin*t) * np.sin(thetas)

	#Unit Position of observer on dust shell in dust shell coords
		xobs = np.sin(thet)*np.cos(phi)
		yobs = np.sin(thet)*np.sin(phi)
		zobs = np.cos(thet)

	#Line of sight velocity anywhere in dust shell (vorb.rhatobs)	
		FracLOS = Vxorb*xobs + Vyorb*yobs + Vzorb*zobs

	#Doppler factor 
		Gams = (1. - bets*bets)**(-0.5)
		Dop = 1./(Gams*(1. - bets*FracLOS))

	# return flux at observer coordinates r, phi, theta (in dust or at earth)
		Fsrc = Lavg/(4.*np.pi*r*r)*(Dop)**(3. - alphnu)



	###-----------------###
	### Compute taudust ###
	###-----------------###
		#xrot = x*np.cos(JJ) + z*np.sin(JJ)
		#zrot = z*np.cos(JJ) - x*np.sin(JJ)
		#throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
		## GET RID OF THIS IF STATEMENT! (did becuase first one catches it)
		
		Qbar=1. ##for now
		tauDust = np.pi*aeff*aeff*Qbar*n0/(1. - p)*(((rofx)/Rd)**(-p) - Rd/rofx) * rofx

		#epsi = 0.0
		#istar=[]
		#Tprof=[]
		#for i in range(len(t)):
		#	istar.append(np.where(Fsrc[i] * np.exp(-tauDust) > RHStable + epsi)[0].max())
		#	Tprof.append(Ttable[istar])

		istar = np.where(Fsrc * np.exp(-tauDust) > RHStable )[0].max()
		Tprof = Ttable[istar]


		#Tprof = (0.25 * Fsrc/sigSB * np.exp(-tauDust) )**(0.25)

	return Tprof










def Fnuint_Shell(ph, thet, nu, t, Dist, Rout, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted form dust
	tem = t - Rd/c*(1. - np.sin(thet)*np.cos(ph))
###----------------------------###
### compute los tau (tauObs) (effective for shell model)   ###
###----------------------------###
	# x = Rd*np.sin(thet)*np.cos(ph)
	# y = Rd*np.sin(thet)*np.sin(ph)
	# z = Rd*np.cos(thet)
	# # doing the integral is faster than lookiup table
	# xe     = Rout*( 1. - (Rd/Rout)*(Rd/Rout) * (  np.cos(thet)*np.cos(thet)  +  np.sin(thet)*np.sin(ph) * np.sin(thet)*np.sin(ph)  )  )**(0.5)
	
	# ##don't integrate if no dust along path
	# if (nDust(xe,y,z, n0, Rd, p, thetT, JJ) == 0.0 and x >= 0.0):
	# 	tauObs = 0.0
	# elif (nDust(xe,y,z, n0, Rd, p, thetT, JJ) == 0.0 and x < 0.0):
	# 	tauObs = np.pi*aeff*aeff * intg.quad(nDust  ,x, 0.0 , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	# else:
	# 	tauObs = np.pi*aeff*aeff * intg.quad(nDust  ,x, xe , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]


	#tauObs = 0.0

	#fint = Qv(nu, nu0, nn) * np.e**(-tauObs) * 2.*h*nu*nu*nu/(c*c)*1./(np.e**(  h*nu/(kb*TDust(tem,Rd, thet, ph, args, RHStable, Ttable))  ) - 1.)	
	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.e**(  h*nu/(kb*TDust(tem,Rd, thet, ph, args, RHStable, Ttable))  ) - 1.)
	fint = fint* Rd*Rd* np.sin(thet) * n0*Rd/(p-1.)


	return np.pi* aeff*aeff/Dist/Dist *fint


def tauObs(x, z, Rout, aeff, n0, Rd, p, thetT, JJ):
	y=0.0
	r    = np.sqrt(x*x + y*y + z*z)
	thet = np.arctan2((x*x + y*y)**(0.5), z)
	ph   = np.arctan2(y, x)

	xe     = Rout*( 1. - (r/Rout)*(r/Rout) * (  np.cos(thet)*np.cos(thet)  +  np.sin(thet)*np.sin(ph) * np.sin(thet)*np.sin(ph)  )  )**(0.5)
	
	#if (type(x) is float):
	return np.pi*aeff*aeff * intg.quad(nDust  ,x, xe , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	#else:
	#	res=[]
	#	for i in range(len(x)):
	#		for j in range(len(z)):
	#				res.append(np.pi*aeff*aeff * intg.quad(nDust  ,x[i], xe[i] , args=(y, z[j], n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0])	

	#	return np.array(res)


def Fnuint_Thick(ph, thet, r, nu, t, Dist, Rout, args, RHStable, Ttable): #, tauGrid):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	# Lavg = args[0]
	# bets = args[1]
	# incl = args[2]
	# Ombin = args[3]
	# alphnu = args[4]    
	# n0 = args[5]
	# Rd = args[6] 
	# p = args[7]
	# thetT = args[8] 
	# JJ = args[9] 
	# aeff = args[10] 
	# nu0 = args[11]
	# nn = args[12]


###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)
## retarded time - time light emitted form dust
	tem =t - r/c*(1. - np.sin(thet)*np.cos(ph))

###----------------------------###
### compute los tau (tauObs)   ###
###----------------------------###
	## doing the integral is faster than lloking it up this way
	xe     = Rout*( 1. - (r/Rout)*(r/Rout) * (  np.cos(thet)*np.cos(thet)  +  np.sin(thet)*np.sin(ph) * np.sin(thet)*np.sin(ph)  )  )**(0.5)
	
	#don't integrate if no dust along path
	if (nDust(xe,y,z, n0, Rd, p, thetT, JJ) == 0.0 and x >= 0.0):
		tauObs = 0.0
	elif (nDust(xe,y,z, n0, Rd, p, thetT, JJ) == 0.0 and x < 0.0):
		tauObs = np.pi*aeff*aeff * intg.quad(nDust  ,x, 0.0 , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	else:
		tauObs = np.pi*aeff*aeff * intg.quad(nDust  ,x, xe , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]


	# #epsi = 2.*Rrout/nn
	# ixl = np.where( x < tauGrid[1])[0].max()
	# ixall = np.where(tauGrid[1][ixl] ==  tauGrid[1])[0]

	# TGz = np.transpose(np.transpose(tauGrid)[ixall])
	# #TGz = np.transpose(TG[ixall])  ## this is all of the z values for given x value
	# it =np.where(z < TGz[2])[0].max()
	# tauObs = tauGrid[0][it]	

	fint = Qv(nu, nu0, nn) * np.e**(-tauObs) * 2.*h*nu*nu*nu/(c*c)*1./(np.e**(  h*nu/(kb*TDust(tem,r, thet, ph, args, RHStable, Ttable))  ) - 1.)
	fint = fint* r*r* np.sin(thet) * nDust(x,y,z, n0, Rd, p, thetT, JJ)


	return np.pi* aeff*aeff/Dist/Dist *fint




# integrate over phi 
def Fnudphi_Shell(thet, nu, t, Dist, Rout, Aargs, RHStable, Ttable):
	#if (type(t) is float):
	return intg.quad(Fnuint_Shell, 0.,2.*np.pi, args=(thet, nu, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	# else:
	# 	res=[]
	# 	i=0
	# 	while (i<len(t)):
	# 		res.append(intg.quad(Fnuint_Shell, 0.,2.*np.pi, args=(thet, nu, t[i], Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0])	
	# 		i += 1
	# 	return np.array(res)

# then int over theta
def Fnu_Shell(nu, t, Dist, Rout, Aargs, RHStable, Ttable):
	#if (type(t) is float):
	return intg.quad(Fnudphi_Shell, 0., np.pi, args=(nu, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]
	#else:
	#	res=[]
	#	i=0
	#	while (i<len(t)):
	#		res.append(intg.quad(Fnudphi_Shell, 0., np.pi, args=(nu, t[i], Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0])	
	#		i += 1
	#	return np.array(res)

# int over freqeuncy
def Fobs_Shell(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable):
	#if (type(t) is float):
	#	return intg.quad(Fnu_Shell, numin, numax, args=(t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]
	#else:
	res=[]
	#i=0
	#while (i<len(t)):
	for i in range(len(t)):
		res.append(intg.quad(Fnu_Shell, numin, numax, args=(t[i], Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo )[0])	
			#i += 1
	return np.array(res)



##FOR MCMC
def magPoint(params, t, THEargs, RHStable, Ttable):
	#beta, cosJJ, Rin, thetT, n0 = params
	cosJJ, Rin, n0 = params
	n0 = n0 * 1.4032428247438431e-09
	Rin = Rin * 2.73213149e+18
	t = t * 86400.
	JJ = np.arccos(cosJJ) ## CAREFUL WITH DOMAIN OF COS
	
	FRel, numin, numax, Dist, Lav, Ombn, alph, pp, Rout,  aeff, nu0, nne, beta, thetT = THEargs
	IncFit = np.arccos(0.07/beta)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
	return -2.5*np.log10(Fobs_Shell(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)







##Now do the same for a thick shell, integrate over r as well...

## integrate over phi 
def Fnudphi_Thick(thet, r, nu, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	return  intg.quad(Fnuint_Thick, 0.,2.*np.pi, args=(thet, r, nu, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	
# then int over theta
def Fnudphidth_Thick(r, nu, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	return intg.quad(Fnudphi_Thick, 0., np.pi, args=(r, nu, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]

# then int over theta
def Fnu_Thick(nu, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	Rd = Aargs[6]
	return intg.quad(Fnudphidth_Thick, Rd, Rout, args=(nu, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]


# int over freqeuncy
def Fobs_Thick(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	if (type(t) is float):
		return intg.quad(Fnu_Thick, numin, numax, args=(t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]
	else:
		res=[]
		i=0
		while (i<len(t)):
			res.append(intg.quad(Fnu_Thick, numin, numax, args=(t[i], Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo )[0])	
			i += 1
		return np.array(res)

