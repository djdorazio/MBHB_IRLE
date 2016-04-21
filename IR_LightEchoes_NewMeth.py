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
#DD CHECKED 4/12/16
def Bv(nu, T):
	return 2.*h*nu*nu*nu/(c*c)*1./(np.exp( h*nu/(kb*T) ) - 1.)
## Dust absorption efficiency
def Qv(nu, nu0, nn):
	# qv = (nu/nu0)**(nn)
	# if (type(qv) is float):
	# 	qv = min(qv, 1.0)
	# 	#if (qv>1.0):
	# 	#	qv=1.0
	# else:
	# 	ii = np.where(qv > 1.0)[0]
	# 	qv[ii] = 1.0
	# return qv
	return np.min( [(nu/nu0)**(nn), 1.*nu/nu])#qv #
	

#DD CHECKED 4/12/16
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



# Torical dust profile DD CHECKED 4/12/16
def nDust(x,y,z, n0, Rd, p, thetT, JJ):
	#nprof = 0.0
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	#rofx  = (xrot*xrot + y*y + zrot*zrot)**(0.5) #same as r
	r  = (x*x + y*y + z*z)**(0.5)
	#r = np.array(r)
	throt = np.arctan2( (xrot*xrot + y*y)**(0.5), zrot)   ##arctan of arg1/arg2 arg1 always positive so btwn 0, pi
	throt = np.array(throt)

	# nprof = 0.0			
	# if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
	# 		nprof = n0*(Rd/r)**(p)
	
	if (type(r) is np.ndarray):
		nprof = 0.0	* r
		for i in range(len(r)):
			if (r[i]>=Rd and throt[i]>=thetT and throt[i]<=(ma.pi - thetT)):
				nprof[i] = n0*(Rd/r[i])**(p)
	else:
		nprof = 0.0			
		if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
			nprof = n0*(Rd/r)**(p)
	

	return nprof



# Torical dust profile
def nD_Schart(x,y,z, n0, Rd, p, thetT, JJ):
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	
	rrot = (xrot*xrot + y*y)**(0.5) # same as r of course
	RzArg = ((rrot*rrot + zrot*zrot)/Rcore/Rcore)**(0.5)
	nprof = n0 * np.e**(  G*Mstr/vt2/Rcore * ( MBH/Mstr/RzArg  + 1./(1.+RzArg)  -  Rcore*(MBH+MstRT)/(2.*Mstr*RT*(1.-p)) *np.pow(r/RT,2.*(p-1.)))  )
	return nprof


## flux from beaming binary as seen by observer at r, theta, phi (centered on binary barycenter)
def Fsrc(t, r, thet, phi, Lavg, bets, incl, Ombin, alphnu):
#Rot by Inc around y axis, Rotation around Bin ang momentum axis by Ombin*t 
#starting point of secindary at t=0
	phis = -ma.pi/2
	thetas = ma.pi/2

	#Make in phase with PG 1302 data (in seconds) # 0.34 from best fit to data
	t = (t - 3600.*24.*1884./(1.+0.2784) * 0.34049274)

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
	return Lavg/(4.*ma.pi*r*r)*(Dop)**(3. - alphnu)


## equation to tabulate RHS and T
#DD CHECKED 4/12/16
def T_RHS(Td, nu0, nn):
	# for for difference in cross sectional area and surface area, pi for isotropic flux from Grain
	RHS = 4.*ma.pi*  (intg.quad(QvBv  ,0., nu0 , args=(Td, nu0, nn) )[0] + intg.quad(Bv  ,nu0 ,numicron*1000., args=(Td) )[0])#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	#return np.log10(RHS)
	return RHS



	


def TDust(t,r,thet,phi,args, RHStable, Ttable):
	#thetT = args[8] 
	#JJ = args[9] 
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	x = r*np.sin(thet)*np.cos(phi)
	y = r*np.sin(thet)*np.sin(phi)
	z = r*np.cos(thet)
	
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	#rofx  = (xrot*xrot + y*y + zrot*zrot)**(0.5) #same as r of course

	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	Tprof = 1.*t/t  ##T=1 is very small

	if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
	###-----------------###
	### COMPUTE Fsrc    ###
	###-----------------###
		#Rot by Inc around y axis, Rotation around Bin ang momentum axis by Ombin*t 
		#starting point of secindary at t=0
		phis = -ma.pi/2.
		thetas = ma.pi/2.


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
		Fsrc = Lavg/(4.*ma.pi*r*r)*(Dop)**(3. - alphnu)



	###-----------------###
	### Compute taudust ###
	###-----------------###
		Qbar=1. ##for now
		tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
		### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
		if (Fsrc * np.exp(-tauDust) > RHStable[len(RHStable)-1] or Fsrc * np.exp(-tauDust) <= RHStable[0]):
			Tprof = 1.
		else:
			istar = np.where(Fsrc * np.exp(-tauDust) <= RHStable )[0].min()
			Tprof = Ttable[istar]

	return Tprof









def TDust_tst(t,r,thet,phi,args, RHStable, Ttable):
	#thetT = args[8] 
	#JJ = args[9] 
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	x = r*np.sin(thet)*np.cos(phi)
	y = r*np.sin(thet)*np.sin(phi)
	z = r*np.cos(thet)
	
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	#rofx  = (xrot*xrot + y*y + zrot*zrot)**(0.5) #same as r of course

	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	Tprof = 1.*t/t  ##T=1 is very small


	if (type(x) is np.ndarray):
		for i in range(len(x)):
			if (r[i]>=Rd and throt[i]>=thetT and throt[i]<=(ma.pi - thetT)):
			###-----------------###
			### COMPUTE Fsrc    ###
			###-----------------###
				#Rot by Inc around y axis, Rotation around Bin ang momentum axis by Ombin*t 
				#starting point of secindary at t=0
				phis = -ma.pi/2.
				thetas = ma.pi/2.


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
				Fsrc = Lavg/(4.*ma.pi*r[i]*r[i])*(Dop)**(3. - alphnu)



			###-----------------###
			### Compute taudust ###
			###-----------------###
				Qbar=1. ##for now
				tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r[i])**(p-1.))
				### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
				if (Fsrc * np.exp(-tauDust) > RHStable[len(RHStable)-1] or Fsrc * np.exp(-tauDust) <= RHStable[0]):
					Tprof = 1.
				else:
					istar = np.where(Fsrc * np.exp(-tauDust) <= RHStable )[0].min()
					Tprof = Ttable[istar]

	else:
		if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
		###-----------------###
		### COMPUTE Fsrc    ###
		###-----------------###
			#Rot by Inc around y axis, Rotation around Bin ang momentum axis by Ombin*t 
			#starting point of secindary at t=0
			phis = -ma.pi/2.
			thetas = ma.pi/2.


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
			Fsrc = Lavg/(4.*ma.pi*r*r)*(Dop)**(3. - alphnu)



		###-----------------###
		### Compute taudust ###
		###-----------------###
			Qbar=1. ##for now
			tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
			### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
			if (Fsrc * np.exp(-tauDust) > RHStable[len(RHStable)-1] or Fsrc * np.exp(-tauDust) <= RHStable[0]):
				Tprof = 1.
			else:
				istar = np.where(Fsrc * np.exp(-tauDust) <= RHStable )[0].min()
				Tprof = Ttable[istar]

	return Tprof







def tauObs(nu, x, y, z, Rout, aeff, n0, Rd, p, thetT, JJ, nu0, nn):

	xe     = (Rout*Rout - (z*z + y*y))**(0.5)
	
	xInt = np.linspace(x, xe, 100.)
	#xInt = np.logspace(x,xe, 100.)
	yInt = nDust(xInt, y, z, n0, Rd, p, thetT, JJ)
	
	#return ma.pi*aeff*aeff * Qv(nu, nu0, nn) * intg.simps(yInt, xInt)
	

	return ma.pi*aeff*aeff * Qv(nu, nu0, nn) * intg.quad(nDust  ,x, xe , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	


## SHELL IS AT r, dust shell starts at Rd
def Fnuint_Shell(ph, thet, nu, r, t, Dist, Rout, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted form dust
	tem = t - r/c*(1. - np.sin(thet)*np.cos(ph))
###----------------------------###
### compute los tau (tauObs) (effective for shell model)   ###
###----------------------------###
	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)

	xe = (Rout*Rout  -  (z*z + y*y))**(0.5)

	#tauobs = ma.pi*aeff*aeff *Qv(nu, nu0, nn)* intg.quad(nDust, x, xe , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	tauobs = tauObs(nu, x, y, z, Rout, aeff, n0, Rd, p, thetT, JJ, nu0, nn)

	#tauObs = 0.0

	#if (  (h*nu/(kb*TDust(tem,Rd, thet, ph, args, RHStable, Ttable))) > 709.7):
	#	fint = 1.654984027680202e+308
	#else:
	#
	#Qbar = 1.
	Surf_nd = 1./(ma.pi * aeff*aeff)#Rde/(p-1.) * ()**(1.-p)
	##RECALL r is the radius of the emititng shell, Rd is the inner edge of the shell
	fint = Qv(nu, nu0, nn) * np.exp(-tauobs) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*TDust(tem,r, thet, ph, args, RHStable, Ttable))  ) - 1.)	
	#fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.e**(  h*nu/(kb*TDust(tem,Rd, thet, ph, args, RHStable, Ttable))  ) - 1.)
	fint = fint* r*r* np.sin(thet) * Surf_nd
	#* n0 * (-(Rd/Rout)**p * Rout + (Rd/r)**p * r)/(-1. + p)
	##^ last term is integral from r to Rout 
	#* n0*Rd/(p-1.) 
	##^last term is the surface density if integrating from rr to Rout for Rout>>Rd

	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist *fint







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
	#xe     = Rout*( 1. - (r/Rout)*(r/Rout) * (  np.cos(thet)*np.cos(thet)  +  np.sin(thet)*np.sin(ph) * np.sin(thet)*np.sin(ph)  )  )**(0.5)
	xe = (Rout*Rout  -  (z*z + y*y))**(0.5)
	#don't integrate if no dust along path
	#if (nDust(xe,y,z, n0, Rd, p, thetT, JJ) == 0.0 and x >= 0.0):
	#	tauObs = 0.0
	#elif (nDust(xe,y,z, n0, Rd, p, thetT, JJ) == 0.0 and x < 0.0):
	#	tauObs = ma.pi*aeff*aeff * intg.quad(nDust  ,x, 0.0 , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	#else:

	#tauObs = ma.pi*aeff*aeff *Qv(nu, nu0, nn)* intg.quad(nDust  ,x, xe , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	tauobs = tauObs(nu, x, y, z, Rout, aeff, n0, Rd, p, thetT, JJ, nu0, nn)

	# #epsi = 2.*Rrout/nn
	# ixl = np.where( x < tauGrid[1])[0].max()
	# ixall = np.where(tauGrid[1][ixl] ==  tauGrid[1])[0]

	# TGz = np.transpose(np.transpose(tauGrid)[ixall])
	# #TGz = np.transpose(TG[ixall])  ## this is all of the z values for given x value
	# it =np.where(z < TGz[2])[0].max()
	# tauObs = tauGrid[0][it]	

	#if (  (h*nu/(kb*TDust(tem,Rd, thet, ph, args, RHStable, Ttable))) > 709.7):
	#	fint = 1.654984027680202e+308
	fint = Qv(nu, nu0, nn) * np.exp(-tauobs) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*TDust(tem,r, thet, ph, args, RHStable, Ttable))  ) - 1.)
	fint = fint* r*r* np.sin(thet) * nDust(x,y,z, n0, Rd, p, thetT, JJ)


	return ma.pi* aeff*aeff/Dist/Dist *fint




# integrate over phi 
def Fnudphi_Shell(thet, nu, r, t, Dist, Rout, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Shell, 0.,2.*ma.pi, args=(thet, nu, r, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	#xInt = np.linspace(0.,2.*ma.pi, 200)
	#yInt = Fnuint_Shell(xInt, thet, nu, r, t, Dist, Rout, Aargs, RHStable, Ttable)
	#return intg.simps(yInt, xInt)
	

	# else:
	# 	res=[]
	# 	i=0
	# 	while (i<len(t)):
	# 		res.append(intg.quad(Fnuint_Shell, 0.,2.*ma.pi, args=(thet, nu, t[i], Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0])	
	# 		i += 1
	# 	return np.array(res)

# then int over theta
def Fnu_Shell(nu, r, t, Dist, Rout, Aargs, RHStable, Ttable):
	#if (type(t) is float):
	return intg.quad(Fnudphi_Shell, 0., ma.pi, args=(nu, r, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]
	#else:
	#	res=[]
	#	i=0
	#	while (i<len(t)):
	#		res.append(intg.quad(Fnudphi_Shell, 0., ma.pi, args=(nu, t[i], Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0])	
	#		i += 1
	#	return np.array(res)

# int over freqeuncy
def Fobs_Shell(numin, numax, r, t, Dist, Rout, Aargs, RHStable, Ttable):
	#if (type(t) is float):
	#	return intg.quad(Fnu_Shell, numin, numax, args=(t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]
	#else:
	res=[]
	#i=0
	#while (i<len(t)):
	for i in range(len(t)):
		res.append(intg.quad(Fnu_Shell, numin, numax, args=(r, t[i], Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo )[0])	
			#i += 1
	return np.array(res)




##FOR MCMC
def magPoint_Shell(params, t, THEargs, RHStable, Ttable, rem_is_Rin):
	if (rem_is_Rin):
		sinJJ, cosTT, Rin, n0 = params
		rem_pls = 0.0
	else:
		sinJJ, cosTT, rem_pls, Rin, n0 = params

	n0 = n0 * 1.4032428247438431e-09
	Rin = Rin * 2.73213149e+18
	rem = Rin + rem_pls * 2.73213149e+18

	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
	
	FRel, numin, numax, Dist, Lav, Ombn, alph, pp, Rout,  aeff, nu0, nne, beta = THEargs
	IncFit = np.arccos(0.067/beta)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
	
	return -2.5*np.log10(Fobs_Shell(numin, numax, rem, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)


	##FOR MCMC
def magPoint_Thick_fmin(params, t, THEargs, RHStable, Ttable):
	#beta, cosJJ, Rin, thetT, n0 = params
	sinJJ, cosTT, Rin, pp, n0 = params
	if (sinJJ<-1.0 or sinJJ >1.0):
		return np.inf
	elif (cosTT<0.0 or cosTT >1.0):
		return np.inf
	else:
		Rin = Rin * 2.73213149e+18
		n0 = n0 * 1.4032428247438431e-09
		t = t * 86400.
		JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
		thetT = np.arccos(cosTT)
		
		#FRel, numin, numax, Dist, Lav, Ombn, alph, Rin, Rout, aeff, nu0, nne, beta = THEargs
		FRel, numin, numax, Dist, Lav, Ombn, alph, Rout, aeff, nu0, nne, beta = THEargs
		IncFit = np.arccos(0.067/beta)

		Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
		#return -2.5*np.log10(Fobs_Shell(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)
		return -2.5*np.log10(Fobs_Thick(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)


def magPoint_Thick(params, t, THEargs, RHStable, Ttable):
	#beta, cosJJ, Rin, thetT, n0 = params
	sinJJ, cosTT, Rin, pp, n0 = params
	Rin = Rin * 2.73213149e+18
	n0 = n0 * 1.4032428247438431e-09
	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
		
	#FRel, numin, numax, Dist, Lav, Ombn, alph, Rin, Rout, aeff, nu0, nne, beta = THEargs
	FRel, numin, numax, Dist, Lav, Ombn, alph, Rout, aeff, nu0, nne, beta = THEargs
	IncFit = np.arccos(0.067/beta)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
	#return -2.5*np.log10(Fobs_Shell(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)
	return -2.5*np.log10(Fobs_Thick(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)




##Now do the same for a thick shell, integrate over r as well...

## integrate over phi 
def Fnudphi_Thick(thet, r, nu, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	return  intg.quad(Fnuint_Thick, 0.,2.*ma.pi, args=(thet, r, nu, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	
# then int over theta
def Fnudphidth_Thick(r, nu, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	return intg.quad(Fnudphi_Thick, 0., ma.pi, args=(r, nu, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]

# then int over theta
def Fnu_Thick(nu, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	Rd = Aargs[6]
	return intg.quad(Fnudphidth_Thick, Rd, Rout, args=(nu, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]

def Fnu_Thick_mult(nu, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	Rd = Aargs[6]
	res=[]
	for i in range (len(nu)):
		res.append(intg.quad(Fnudphidth_Thick, Rd, Rout, args=(nu[i], t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0])
	return np.array(res)




# int over freqeuncy
def Fobs_Thick(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	#if (type(t) is float):
	#	return intg.quad(Fnu_Thick, numin, numax, args=(t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]
	#else:
	res=[]
	#i=0
	for i in range (len(t)):
		res.append(intg.quad(Fnu_Thick, numin, numax, args=(t[i], Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo )[0])	
		#i += 1
	return np.array(res)

def Fobs_Thick_n0(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	res=[]
	for i in range (len(Aargs)):
		res.append(intg.quad(Fnu_Thick, numin, numax, args=(t, Dist, Rout, Aargs[i], RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo )[0])	
		#i += 1
	return np.array(res)

