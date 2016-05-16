import numpy as np
import math as ma
#import numexpr as ne
import scipy as sc

#from scipy.optimize import brentq #fmin

import scipy.integrate as intg
import scipy.signal as sgn
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
	qv = (nu/nu0)**(nn)
#	if (type(qv) is float):
	qv = min(qv, 1.0)
#	else:
#		ii = np.where(qv > 1.0)[0]
#		qv[ii] = 1.0
	return qv
	#return np.min( [(nu/nu0)**(nn), 1.*nu/nu])
#Qv = np.vectorize(Qv)

#DD CHECKED 4/12/16
def QvBv(nu, T, nu0, nn):
	qv = (nu/nu0)**(nn)
	#if (type(qv) is float):
	qv = min(qv, 1.0)
		#if (qv>1.0):
		#	qv=1.0
	# else:
	# 	ii = np.where(qv > 1.0)[0]
	# 	qv[ii] = 1.0
	return 2.*h*nu*nu*nu/(c*c)*1./(np.exp(h*nu/(kb*T)) - 1.) * qv
#QvBv = np.vectorize(QvBv)


# Torical dust profile
def nDust(x,y,z, n0, Rd, p, thetT, JJ):
	#nprof = 0.0
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	#rofx  = (xrot*xrot + y*y + zrot*zrot)**(0.5) #same as r
	r  = (x*x + y*y + z*z)**(0.5)

	throt = np.arctan2( (xrot*xrot + y*y)**(0.5), zrot)   ##arctan of arg1/arg2 arg1 always positive so btwn 0, pi
	throt = np.array(throt)

	
	# if (type(r) is np.ndarray):
	# 	nprof = 0.0	* r
	# 	for i in range(len(r)):
	# 		if (r[i]>=Rd and throt[i]>=thetT and throt[i]<=(ma.pi - thetT)):
	# 			nprof[i] = n0*(Rd/r[i])**(p)
	# else:
	nprof = 0.0			
	if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
		nprof = n0*(Rd/r)**(p)
	

	return nprof

#nDust = np.vectorize(nDust, excluded = (3,4,5,6,7), cache=True)




# Torical dust profile
def nprof(r, n0, Rd, p):
	return n0*(Rd/r)**(p)

def nDust_pcwse(x,y,z, n0, Rd, p, thetT, JJ):
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	r  = (x*x + y*y + z*z)**(0.5)

	throt = np.arctan2( (xrot*xrot + y*y)**(0.5), zrot)   ##arctan of arg1/arg2 arg1 always positive so btwn 0, pi


	#nprof = n0*(Rd/r)**(p)
	#BoxCar = Heaviside[throt - thetT] - Heaviside[throt - (ma.pi - thetT)]
	#return nprof * Heaviside[r-Rd] * BoxCar

	#nd = np.piecewise(r, [r < Rd, r >= Rd], [lambda r:0.0, lambda r:n0*(Rd/r)**(p)])
	nd = np.piecewise(throt, [throt>=(np.pi - thetT), throt<=(np.pi - thetT)], [lambda throt:0.0, lambda throt:n0])
	nd = np.piecewise(throt, [throt<thetT, throt>=thetT], [lambda throt:0.0, lambda throt:nd])

#or throt>=(np.pi - thetT)
#and throt<=(np.pi - thetT)
	return nd






## equation to tabulate RHS and T
def T_RHS(Td, nu0, nn):
	# for for difference in cross sectional area and surface area, pi for isotropic flux from Grain
	RHS = 4.*ma.pi*  (intg.quad(QvBv  ,0., nu0 , args=(Td, nu0, nn) )[0] + intg.quad(Bv  ,nu0 ,numicron*1000., args=(Td) )[0])#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	#return np.log10(RHS)
	return RHS













def tauObs(nu, x, y, z, Rout, aeff, n0, Rd, p, thetT, JJ, nu0, nn):
	xe     = (Rout*Rout - (z*z + y*y))**(0.5)
## SIMPSON
	#xInt = np.linspace(x, xe, 100.)
	#xInt = np.logspace(x,xe, 100.)
	#yInt = nDust(xInt, y, z, n0, Rd, p, thetT, JJ)
	#return ma.pi*aeff*aeff * Qv(nu, nu0, nn) * intg.simps(yInt, xInt)	

## GAUSSIAN QUAD
	return ma.pi*aeff*aeff * Qv(nu, nu0, nn) * intg.quad(nDust  ,x, xe , args=(y, z, n0, Rd, p, thetT, JJ) , epsabs=myabs, epsrel=myrel, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	


def tauObs_Shell(nu, x, y, z, aeff, n0, Rd, p, thetT, JJ, nu0, nn):
	Rout = Rd
	##xe is the x value on teh near side of the dust torus
	#xe     = (Rout*Rout - (z*z + y*y))**(0.5)
	xe = (x*x)**0.5 

	#xin    = (Rd*Rd - (z*z + y*y))**(0.5)
	#if (x > xe):
	#	print "ERROR in xe"
	#	import sys
	#	sys.exit(0)

##INNER SHELL
	# if ndust is zero on the near side of the dust, then tau=0
	if (nDust(xe,y,z, n0, Rd, p, thetT, JJ) == 0.0):#(z*z + y*y) <= Rd*Rd):
		tau = 0.00000000#nDust(xe,y,z, n0, Rd, p, thetT, JJ)
	else:
		tau = 100000000.0
	
	return tau
















####################################################
###ISO CASES    
####################################################
def Fsrc_Iso(t, r, Lavg, Amp, Ombin, t0):
	# Amp is in fraction of total flux, t and t0 in seconds
	#return Lavg/(4.*ma.pi*r*r)* ( 1. + Amp*np.sin(Ombin*(t - t0))  )
	return Lavg/(4.*ma.pi*r*r)* ( 1. + Amp*np.sin(Ombin*t - t0)  )  #make t0 pahse not time






def TDust_Iso(t,r,thet, ph, args, RHStable, Ttable):
	Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)
	
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	#rofx  = (xrot*xrot + y*y + zrot*zrot)**(0.5) #same as r of course

	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	Td = 1.*r/r  ##T=1 is very small

	if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
		###-----------------###
		### COMPUTE Fsrc    ###
		###-----------------###

		Fsrc = Fsrc_Iso(t, r, Lavg, Amp, Ombin, t0)


		###-----------------###
		### Compute taudust ###
		###-----------------###
		Qbar=1. 
		tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
		### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
		LHS = Fsrc * np.exp(-tauDust)
		# if (type(Fsrc) is np.ndarray):
		# 	if (len(Fsrc) > len(RHStable)):
		# 		LHS = sgn.resample(LHS, len(RHStable))
		# 	elif (len(Fsrc) < len(RHStable)):
		# 		RHStable = sgn.resample(RHStable, len(Fsrc))
		RHS_mx = RHStable[len(RHStable)-1]
		RHS_mn = RHStable[0]

		if (LHS > RHS_mx or LHS <= RHS_mn):
			Td = 1.
		else:
			istar = np.where( LHS <= RHStable )[0].min()
			Td = Ttable[istar]

	return Td
	


## SHELL TORUS OPTHIN Fnu for doppler beaming case, optically thick shell torus
def Fnuint_Ring_Iso(ph, nu, t, Dist, args, RHStable, Ttable):
	Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args

###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## Coordinates of the dust Ring
	yy = np.sign(np.cos(ph)) * np.cos(JJ) * (1. +  (np.tan(ph)/np.cos(JJ))**2 )**(0.5) 
	xx = -np.sin(JJ)
	Th_ring = np.arctan2( yy, xx ) - np.pi/2. * (np.sign(np.cos(ph)) - 1.)
	#Th_ring = ma.pi/2.
## retarded time - time light emitted form dust
	tem = t - Rd/c*( 1. - np.sin(Th_ring)*np.cos(ph) ) 

	
	# Tdust for ISO source
	args[7] = 0.0 #set thea_T = 0 for Ring
	Tdust = TDust_Iso(tem, Rd, Th_ring, ph, args, RHStable, Ttable)
	# surface density in optically thick limit
	l_d = 1./(2.*aeff)
	# Rd is the inner edge of the shell
	fint =  2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = fint* Rd*Rd* np.sin(Th_ring) * l_d #*Qv(nu, nu0, nn) * 
	

	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist *fint


### This is the spherical shell case for theta_T = 0
## SHELL TORUS OPTHIN Fnu for doppler beaming case, optically thick shell torus
def Fnuint_Shell_OptThin_Iso(ph, thet, nu, t, Dist, args, RHStable, Ttable):
	Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args

###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted form dust
	tem = t - Rd/c*(1. - np.sin(thet)*np.cos(ph))

	# Tdust for doppler source
	Tdust = TDust_Iso(tem, Rd, thet, ph,  args, RHStable, Ttable)
	# surface density in optically thick limit
	Surf_nd = 1./(ma.pi * aeff*aeff)
	# Rd is the inner edge of the shell
	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = fint* Rd*Rd* np.sin(thet) * Surf_nd
	

	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist *fint


## Fnu for doppler beaming case, optically thick shell torus
def Fnuint_Shell_OptThick_Iso(ph, thet, nu, t, Dist, args, RHStable, Ttable):
	Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args

###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted form dust
	tem = t - Rd/c*(1. - np.sin(thet)*np.cos(ph))
###----------------------------###
### compute los tau (tauObs) (effective for shell model)   ###
###----------------------------###
	x = Rd*np.sin(thet)*np.cos(ph)
	y = Rd*np.sin(thet)*np.sin(ph)
	z = Rd*np.cos(thet)

	#In optically thick case blocck light which hits interverning dust
	tauobs = tauObs_Shell(nu, x, y, z, aeff, n0, Rd, p, thetT, JJ, nu0, nn)

	# Tdust for doppler source
	Tdust = TDust_Iso(tem, Rd, thet, ph, args, RHStable, Ttable)
	# surface density in optically thick limit
	Surf_nd = 1./(ma.pi * aeff*aeff)
	# Rd is the inner edge of the shell
	fint = Qv(nu, nu0, nn) * np.exp(-tauobs) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = fint* Rd*Rd* np.sin(thet) * Surf_nd
	

	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist *fint



def Fnuint_Thick_Iso(ph, thet, r, nu, t, Dist, args, RHStable, Ttable): #, tauGrid):
	Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	## where UV penetrates out to (tau_UV=1)
	Rout = Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) + 0.01*Rd



###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)
## retarded time - time light emitted form dust
	tem = t - r/c*(1. - np.sin(thet)*np.cos(ph))



	# Tdust for doppler source
	Tdust = TDust_Iso(tem, r, thet, ph, args, RHStable, Ttable)

	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)
	fint = fint* r*r* np.sin(thet) * nDust(x,y,z, n0, Rd, p, thetT, JJ)


	return ma.pi* aeff*aeff/Dist/Dist *fint


####################################################
### END ISO CASES    
####################################################
























####################################################
###Doppler CASES    
####################################################
## flux from beaming binary as seen by observer at r, theta, phi (centered on binary barycenter)
def Fsrc_Dop_PG(t, r, thet, phi, Lavg, bets, incl, Ombin, alphnu):
#Rot by Inc around y axis, Rotation around Bin ang momentum axis by Ombin*t 
#starting point of secindary at t=0

	## best fit to PG 1302 optical
	#[fraction of lum in optical, beta, phase, inclination in radians, obs period in days]
	Fsrc_p_opt = [ 5.98144879e-02,   6.12468791e-02,   6.55067929e-01,  -3.28334799e-04, 1.87091995e+03]

	Ombin = 2.*ma.pi/(Fsrc_p_opt[4]*24.*3600.)* (1.+0.2784)
	t = t-Fsrc_p_opt[2]*2.*ma.pi/Ombin

	bets = Fsrc_p_opt[1]
	incl = Fsrc_p_opt[2]


	phis = -ma.pi/2
	thetas = ma.pi/2

	#Make in phase with PG 1302 data (in seconds) # 0.34 from best fit to data
	#t = (t - 3600.*24.*1884./(1.+0.2784) * 0.34049274)

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



def Fsrc_Dop(t, r, thet, phi, Lavg, bets, incl, Ombin, alphnu):
#Rot by Inc around y axis, Rotation around Bin ang momentum axis by Ombin*t 
#starting point of secindary at t=0
	phis = -ma.pi/2
	thetas = ma.pi/2

	#Make in phase with PG 1302 data (in seconds) # 0.34 from best fit to data
	#t = (t - 3600.*24.*1884./(1.+0.2784) * 0.34049274)

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




def TDust_Dop_pcw(t,r,thet,phi,args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	x = r*np.sin(thet)*np.cos(phi)
	y = r*np.sin(thet)*np.sin(phi)
	z = r*np.cos(thet)
	
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	#rofx  = (xrot*xrot + y*y + zrot*zrot)**(0.5) #same as r of course

	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	

	#Tprof = 1.*t/t  ##T=1 is very small

	#if (r>=Rd and throt>=thetT and throt<=(ma.pi - thetT)):
	###-----------------###
	### COMPUTE Fsrc    ###
	###-----------------###

	Fsrc = Fsrc_Dop(t, r, thet, phi, Lavg, bets, incl, Ombin, alphnu)


	###-----------------###
	### Compute taudust ###
	###-----------------###
	Qbar=1. ##for now
	tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
	### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
	
	LHS = Fsrc * np.exp(-tauDust)
	if (type(Fsrc) is np.ndarray):
		if (len(Fsrc) > len(RHStable)):
			LHS = sgn.resample(LHS, len(RHStable))
		elif (len(Fsrc) < len(RHStable)):
			RHStable = sgn.resample(RHStable, len(Fsrc))

	#else:


	RHS_mx = RHStable[len(RHStable)-1]
	RHS_mn = RHStable[0]

	#Td = np.piecewise(r, [LHS > RHS_mx, LHS <= RHS_mn], [lambda LHS:1.0, lambda LHS:Ttable[np.where(LHS <= RHStable )[0].min()]])

	#if (LHS > RHS_mx, LHS <= RHS_mn):
	#	Tprof = 1.
	#else:
	#for i in range():
	#LHS = max(LHS, RHS_mn)
	#LHS = min(LHS, RHS_mn)

	#iNoT = np.where(LHS > RHS_mx or LHS <= RHS_mn)


	istar = np.where( LHS <= RHStable )[0].min()
	Td = Ttable[istar]



	#Td = np.piecewise(r, [r < Rd, r >= Rd], [lambda r:0.0, lambda r:Td])
	#Td = np.piecewise(throt, [throt<thetT or throt>=(np.pi - thetT), throt>=thetT and throt<=(np.pi - thetT)], [lambda throt:0.0, lambda throt:Td])


	return Td




def TDust_Dop(t,r,thet,phi, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	x = r*np.sin(thet)*np.cos(phi)
	y = r*np.sin(thet)*np.sin(phi)
	z = r*np.cos(thet)
	
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	

	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	

	Td = 1.*r/r ##T=1 is very small
	# if (type(r) is np.ndarray):
	# 	for i in range(len(r)):
	# 		if ( r[i]>=Rd and  throt[i]>=thetT and  throt[i]<=(ma.pi - thetT)):
	# 			###-----------------###
	# 			### COMPUTE Fsrc    ###
	# 			###-----------------###
	# 			Fsrc = Fsrc_Dop(t, r[i], thet[i], phi[i], Lavg, bets, incl, Ombin, alphnu)
	# 			###-----------------###
	# 			### Compute taudust ###
	# 			###-----------------###
	# 			Qbar=1. ##for now
	# 			tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r[i])**(p-1.))
	# 			### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
	# 			LHS = Fsrc * np.exp(-tauDust)
	# 			if (type(Fsrc) is np.ndarray):
	# 				if (len(Fsrc) > len(RHStable)):
	# 					LHS = sgn.resample(LHS, len(RHStable))
	# 				elif (len(Fsrc) < len(RHStable)):
	# 					RHStable = sgn.resample(RHStable, len(Fsrc))
	# 			RHS_mx = RHStable[len(RHStable)-1]
	# 			RHS_mn = RHStable[0]

	# 			#Td = np.piecewise(r, [LHS > RHS_mx, LHS <= RHS_mn], [lambda LHS:1.0, lambda LHS:Ttable[np.where(LHS <= RHStable )[0].min()]])
	# 			#if (LHS > RHS_mx, LHS <= RHS_mn):
	# 			#	Td = 1.
	# 			#else:
	# 			istar = np.where( LHS <= RHStable )[0].min()
	# 			Td[i] = Ttable[istar]
	# else:
	if ( r>=Rd and  throt>=thetT and  throt<=(ma.pi - thetT)):
		###-----------------###
		### COMPUTE Fsrc    ###
		###-----------------###
		Fsrc = Fsrc_Dop(t, r, thet, phi, Lavg, bets, incl, Ombin, alphnu)
		###-----------------###
		### Compute taudust ###
		###-----------------###
		Qbar=1. ##for now
		tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
		### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
		LHS = Fsrc * np.exp(-tauDust)
		# if (type(Fsrc) is np.ndarray):
		# 	if (len(Fsrc) > len(RHStable)):
		# 		LHS = sgn.resample(LHS, len(RHStable))
		# 	elif (len(Fsrc) < len(RHStable)):
		# 		RHStable = sgn.resample(RHStable, len(Fsrc))
		RHS_mx = RHStable[len(RHStable)-1]
		RHS_mn = RHStable[0]

		#Td = np.piecewise(r, [LHS > RHS_mx, LHS <= RHS_mn], [lambda LHS:1.0, lambda LHS:Ttable[np.where(LHS <= RHStable )[0].min()]])
		if (LHS > RHS_mx or LHS <= RHS_mn):
			Td = 1.
		else:
			istar = np.where( LHS <= RHStable )[0].min()
			Td = Ttable[istar]

	return Td




def TDust_Dop_PG(t,r,thet,phi, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	x = r*np.sin(thet)*np.cos(phi)
	y = r*np.sin(thet)*np.sin(phi)
	z = r*np.cos(thet)
	
	xrot = x*np.cos(JJ) + z*np.sin(JJ)
	zrot = z*np.cos(JJ) - x*np.sin(JJ)
	

	throt = np.arctan2((xrot*xrot + y*y)**(0.5), zrot)
	

	Td = 1.*r/r ##T=1 is very small
	# if (type(r) is np.ndarray):
	# 	for i in range(len(r)):
	# 		if ( r[i]>=Rd and  throt[i]>=thetT and  throt[i]<=(ma.pi - thetT)):
	# 			###-----------------###
	# 			### COMPUTE Fsrc    ###
	# 			###-----------------###
	# 			Fsrc = Fsrc_Dop(t, r[i], thet[i], phi[i], Lavg, bets, incl, Ombin, alphnu)
	# 			###-----------------###
	# 			### Compute taudust ###
	# 			###-----------------###
	# 			Qbar=1. ##for now
	# 			tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r[i])**(p-1.))
	# 			### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
	# 			LHS = Fsrc * np.exp(-tauDust)
	# 			if (type(Fsrc) is np.ndarray):
	# 				if (len(Fsrc) > len(RHStable)):
	# 					LHS = sgn.resample(LHS, len(RHStable))
	# 				elif (len(Fsrc) < len(RHStable)):
	# 					RHStable = sgn.resample(RHStable, len(Fsrc))
	# 			RHS_mx = RHStable[len(RHStable)-1]
	# 			RHS_mn = RHStable[0]

	# 			#Td = np.piecewise(r, [LHS > RHS_mx, LHS <= RHS_mn], [lambda LHS:1.0, lambda LHS:Ttable[np.where(LHS <= RHStable )[0].min()]])
	# 			#if (LHS > RHS_mx, LHS <= RHS_mn):
	# 			#	Td = 1.
	# 			#else:
	# 			istar = np.where( LHS <= RHStable )[0].min()
	# 			Td[i] = Ttable[istar]
	# else:
	if ( r>=Rd and  throt>=thetT and  throt<=(ma.pi - thetT)):
		###-----------------###
		### COMPUTE Fsrc    ###
		###-----------------###
		Fsrc = Fsrc_Dop_PG(t, r, thet, phi, Lavg, bets, incl, Ombin, alphnu)
		###-----------------###
		### Compute taudust ###
		###-----------------###
		Qbar=1. ##for now
		tauDust = ma.pi*aeff*aeff*Qbar*n0/(p-1.)*  Rd *( 1 -  (Rd/r)**(p-1.))
		### if flux is greater than RHS max at which T > Tsub~2000K, then dust sublimates
		LHS = Fsrc * np.exp(-tauDust)
		# if (type(Fsrc) is np.ndarray):
		# 	if (len(Fsrc) > len(RHStable)):
		# 		LHS = sgn.resample(LHS, len(RHStable))
		# 	elif (len(Fsrc) < len(RHStable)):
		# 		RHStable = sgn.resample(RHStable, len(Fsrc))
		RHS_mx = RHStable[len(RHStable)-1]
		RHS_mn = RHStable[0]

		#Td = np.piecewise(r, [LHS > RHS_mx, LHS <= RHS_mn], [lambda LHS:1.0, lambda LHS:Ttable[np.where(LHS <= RHStable )[0].min()]])
		if (LHS > RHS_mx or LHS <= RHS_mn):
			Td = 1.
		else:
			istar = np.where( LHS <= RHStable )[0].min()
			Td = Ttable[istar]

	return Td	


#TDust_Dop = np.vectorize(TDust_Dop, excluded=(4,5,6), cache=True)#excluded=("args", "RHStable", "Ttable"))


def Fnuint_Ring_Dop(ph, nu, t, Dist, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args

###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## Coordinates of the dust Ring
	yy = np.sign(np.cos(ph)) * np.cos(JJ) * (1. +  (np.tan(ph)/np.cos(JJ))**2 )**(0.5) 
	xx = -np.sin(JJ)
	Th_ring = np.arctan2( yy, xx ) - np.pi/2. * (np.sign(np.cos(ph)) - 1.)
	#Th_ring = ma.pi/2.
## retarded time - time light emitted form dust
	tem = t - Rd/c*( 1. - np.sin(Th_ring)*np.cos(ph) ) 

	
	# Tdust for ISO source
	args[8] = 0.0 #set thea_T = 0 for Ring/Sphere
	Tdust = TDust_Dop(tem, Rd, Th_ring, ph, args, RHStable, Ttable)
	# surface density in optically thick limit
	l_d = 1./(2.*aeff)
	# Rd is the inner edge of the shell
	fint =  2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = fint* Rd*Rd* np.sin(Th_ring) * l_d #*Qv(nu, nu0, nn) * 
	

	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist *fint

## Fnu for doppler beaming case, optically thick shell torus
def Fnuint_Shell_OptThin_Dop(ph, thet, nu, t, Dist, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args

###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted form dust
	tem = t - Rd/c*(1. - np.sin(thet)*np.cos(ph))
###----------------------------###
### compute los tau (tauObs) (effective for shell model)   ###
###----------------------------###
	#x = Rd*np.sin(thet)*np.cos(ph)
	#y = Rd*np.sin(thet)*np.sin(ph)
	#z = Rd*np.cos(thet)

	
	# Tdust for doppler source
	Tdust = TDust_Dop(tem,Rd, thet, ph, args, RHStable, Ttable)
	# surface density in optically thick limit
	Surf_nd = 1./(ma.pi * aeff*aeff)
	# Rd is the inner edge of the shell
	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = fint* Rd*Rd* np.sin(thet) * Surf_nd
	

	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist *fint


def Fnuint_Shell_OptThin_Dop_PG(ph, thet, nu, t, Dist, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args

###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted form dust
	tem = t - Rd/c*(1. - np.sin(thet)*np.cos(ph))
###----------------------------###
### compute los tau (tauObs) (effective for shell model)   ###
###----------------------------###
	#x = Rd*np.sin(thet)*np.cos(ph)
	#y = Rd*np.sin(thet)*np.sin(ph)
	#z = Rd*np.cos(thet)

	
	# Tdust for doppler source
	Tdust = TDust_Dop_PG(tem,Rd, thet, ph, args, RHStable, Ttable)
	# surface density in optically thick limit
	Surf_nd = 1./(ma.pi * aeff*aeff)
	# Rd is the inner edge of the shell
	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = fint* Rd*Rd* np.sin(thet) * Surf_nd
	

	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist *fint


## Fnu for doppler beaming case, optically thick shell torus
def Fnuint_Shell_OptThick_Dop(ph, thet, nu, t, Dist, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args

###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted form dust
	tem = t - Rd/c*(1. - np.sin(thet)*np.cos(ph))
###----------------------------###
### compute los tau (tauObs) (effective for shell model)   ###
###----------------------------###
	x = Rd*np.sin(thet)*np.cos(ph)
	y = Rd*np.sin(thet)*np.sin(ph)
	z = Rd*np.cos(thet)

	#In optically thick case blocck light which hits interverning dust
	#tauobs = tauObs_Shell(nu, x, y, z, aeff, n0, Rd, p, thetT, JJ, nu0, nn)

	# Tdust for doppler source
	Tdust = TDust_Dop(tem, Rd, thet, ph, args, RHStable, Ttable)
	# surface density in optically thick limit
	Surf_nd = 1./(ma.pi * aeff*aeff)
	# Rd is the inner edge of the shell
	#fint = Qv(nu, nu0, nn) * np.exp(-tauobs) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = fint* Rd*Rd* np.sin(thet) * Surf_nd
	

	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist *fint



def Fnuint_Shell_OptThick_Dop_PG(ph, thet, nu, t, Dist, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args

###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
## retarded time - time light emitted form dust
	tem = t - Rd/c*(1. - np.sin(thet)*np.cos(ph))
###----------------------------###
### compute los tau (tauObs) (effective for shell model)   ###
###----------------------------###
	x = Rd*np.sin(thet)*np.cos(ph)
	y = Rd*np.sin(thet)*np.sin(ph)
	z = Rd*np.cos(thet)

	#In optically thick case blocck light which hits interverning dust
	#tauobs = tauObs_Shell(nu, x, y, z, aeff, n0, Rd, p, thetT, JJ, nu0, nn)

	# Tdust for doppler source
	Tdust = TDust_Dop_PG(tem, Rd, thet, ph, args, RHStable, Ttable)
	# surface density in optically thick limit
	Surf_nd = 1./(ma.pi * aeff*aeff)
	# Rd is the inner edge of the shell
	#fint = Qv(nu, nu0, nn) * np.exp(-tauobs) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)	
	fint = fint* Rd*Rd* np.sin(thet) * Surf_nd
	

	# pi for uniform emitting dust grain
	return ma.pi* aeff*aeff/Dist/Dist *fint





def Fnuint_Thick_Dop(ph, thet, r, nu, t, Dist, args, RHStable, Ttable): #, tauGrid):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	## where UV penetrates out to (tau_UV=1)
	Rout = Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) + 0.01*Rd



###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)
## retarded time - time light emitted form dust
	tem = t - r/c*(1. - np.sin(thet)*np.cos(ph))



	# Tdust for doppler source
	Tdust = TDust_Dop(tem, r, thet, ph, args, RHStable, Ttable)

	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)
	fint = fint* r*r* np.sin(thet) * nDust(x,y,z, n0, Rd, p, thetT, JJ)
	#fint = fint * nDust(x,y,z, n0, Rd, p, thetT, JJ)


	return ma.pi* aeff*aeff/Dist/Dist *fint




def Fnuint_Thick_Dop_PG(ph, thet, r, nu, t, Dist, args, RHStable, Ttable): #, tauGrid):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	## where UV penetrates out to (tau_UV=1)
	Rout = Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) + 0.01*Rd



###----------------------------###
### SETUP COORDS TO INTEGRATE  ###
###----------------------------###
	x = r*np.sin(thet)*np.cos(ph)
	y = r*np.sin(thet)*np.sin(ph)
	z = r*np.cos(thet)
## retarded time - time light emitted form dust
	tem = t - r/c*(1. - np.sin(thet)*np.cos(ph))



	# Tdust for doppler source
	Tdust = TDust_Dop_PG(tem, r, thet, ph, args, RHStable, Ttable)

	fint = Qv(nu, nu0, nn) * 2.*h*nu*nu*nu/(c*c)*1./(np.exp(  h*nu/(kb*Tdust)  ) - 1.)
	fint = fint* r*r* np.sin(thet) * nDust(x,y,z, n0, Rd, p, thetT, JJ)
	#fint = fint * nDust(x,y,z, n0, Rd, p, thetT, JJ)


	return ma.pi* aeff*aeff/Dist/Dist *fint




####################################################
###   END Doppler CASES    
####################################################
















####################################################
###   INTEGRATION   
####################################################


####################################################
###   MC INTEGRATION   
####################################################


def Thick_Dop_MCInt(N, numin, numax, t, Dist, args, RHStable, Ttable):
	Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	## where UV penetrates out to (tau_UV=1)
	Rtau1 = Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) + 0.01*Rd
 
	#r_rnd  = Rd + np.random.rand(N)* (Rtau1-Rd)
	#Rtau1=1
	#Rd=0.0
	r_rnd  = ( (Rtau1*Rtau1*Rtau1 - Rd*Rd*Rd)*np.random.rand(N) +  Rd*Rd*Rd)**(1./3.)

	th_rnd = np.arccos(1. -  2.*np.random.rand(N))
	ph_rnd = np.random.rand(N) * 2.*np.pi
	nu_rnd = numin + np.random.rand(N) * (numax - numin)

	#NTot = N
	Vol = 4.*ma.pi/3. * (Rtau1*Rtau1*Rtau1 - Rd*Rd*Rd) * (numax-numin)  #4.*ma.pi * (Rtau1 - Rd)#
	#Vol = 4.*ma.pi/3. * 
	#return Vol/N * np.sum(r_rnd*r_rnd)
	return Vol/N * np.sum( Fnuint_Thick_Dop(ph_rnd, th_rnd, r_rnd, nu_rnd, t, Dist, args, RHStable, Ttable) ) 



####################################################
###   GAUSSIAN QUADRATURE
####################################################
##################
### F_ISO
##################
#########
###Sphere
#########
def FThnu_Sphere_Iso_QuadInt(thet, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Shell_OptThin_Iso, 0.,2.*ma.pi, args=(thet, nu, t, Dist, Aargs, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,full_output=fo  )[0]

def Fnu_Sphere_Iso_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThnu_Sphere_Iso_QuadInt, 0., ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]

def F_Sphere_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	Aargs[7] = 0.0
	return intg.quad(Fnu_Sphere_Iso_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]



#########
###RING
#########
#Fnuint_Ring_Iso(ph, thet, nu, t, Dist, args, RHStable, Ttable)
def Fnu_Ring_Iso_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Ring_Iso, 0.,2.*ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable),  full_output=fo  )[0]

def F_Ring_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnu_Ring_Iso_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable),  full_output=fo  )[0]
#F_Ring_Iso_QuadInt = np.vectorize(F_Ring_Iso_QuadInt, excluded=(0,1,3,4,5,6), cache=True)#excluded=("args", "RHStable", "Ttable"))




#########
###Torus Shell - OPT THIN
#########
def F_ShTorOptThin_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnu_Sphere_Iso_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]




#########
###Torus Shell - OPT THICK
#########
def FThnu_ShTorOptThick_Iso_QuadInt(thet, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Shell_OptThick_Iso, 0.,2.*ma.pi, args=(thet, nu, t, Dist, Aargs, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,full_output=fo  )[0]

def Fnu_ShTorOptThick_Iso_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThnu_ShTorOptThick_Iso_QuadInt, 0., ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]

def F_ShTorOptThick_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnu_ShTorOptThick_Iso_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]



#########
###THICK
#########
def FThRnu_Thick_Iso_QuadInt(thet, r, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Thick_Iso, 0.,2.*ma.pi, args=(thet, r, nu, t, Dist, args, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def FRnu_Thick_Iso_QuadInt(r, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThRnu_Thick_Iso_QuadInt, 0., ma.pi, args=(r, nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def Fnu_Thick_Iso_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	#Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	n0 = Aargs[4]
	Rd = Aargs[5]
	p  = Aargs[6]
	aeff = Aargs[9]
	Rtau1 = Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) + 0.01*Rd
	return intg.quad(FRnu_Thick_Iso_QuadInt, Rd, Rtau1, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def F_Thick_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnu_Thick_Iso_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]









##################
### F_DOPPPLER
##################
#########
###Sphere
#########
def FThnu_Sphere_Dop_QuadInt(thet, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Shell_OptThin_Dop, 0.,2.*ma.pi, args=(thet, nu, t, Dist, Aargs, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,full_output=fo  )[0]
		
def Fnu_Sphere_Dop_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThnu_Sphere_Dop_QuadInt, 0., ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]

def F_Sphere_Dop_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	Aargs[8] = 0.0 #set thea_T = 0 for Ring/Sphere
	return intg.quad(Fnu_Sphere_Dop_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]




#########
###RING
#########
#Fnuint_Ring_Iso(ph, thet, nu, t, Dist, args, RHStable, Ttable)
def Fnu_Ring_Dop_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Ring_Dop, 0.,2.*ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable),  full_output=fo  )[0]

def F_Ring_Dop_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnu_Ring_Dop_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable),  full_output=fo  )[0]
#F_Ring_Iso_QuadInt = np.vectorize(F_Ring_Iso_QuadInt, excluded=(0,1,3,4,5,6), cache=True)#excluded=("args", "RHStable", "Ttable"))


#########
###Shell Torus Opt Thin
#########
def F_ShTorOptThin_Dop_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnu_Sphere_Dop_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]


#########
###Shell Torus Opt Thick
#########
def FThnu_ShTorOptThick_Dop_QuadInt(thet, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Shell_OptThick_Dop, 0.,2.*ma.pi, args=(thet, nu, t, Dist, Aargs, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,full_output=fo  )[0]
		
def Fnu_ShTorOptThick_Dop_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThnu_ShTorOptThick_Dop_QuadInt, 0., ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]

def F_ShTorOptThick_Dop_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	#return intg.quad(Fnu_ShTorOptThick_Dop_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]
	 res=[]
	 for i in range(len(t)):
	 	res.append(intg.quad(Fnu_ShTorOptThick_Dop_QuadInt, numin, numax, args=(t[i], Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0] )	
	 return np.array(res)

###
def FThnu_ShTorOptThick_Dop_QuadInt_PG(thet, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Shell_OptThick_Dop_PG, 0.,2.*ma.pi, args=(thet, nu, t, Dist, Aargs, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,full_output=fo  )[0]
		
def Fnu_ShTorOptThick_Dop_QuadInt_PG(nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThnu_ShTorOptThick_Dop_QuadInt_PG, 0., ma.pi, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]

def F_ShTorOptThick_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	#return intg.quad(Fnu_ShTorOptThick_Dop_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0]
	 res=[]
	 for i in range(len(t)):
	 	res.append(intg.quad(Fnu_ShTorOptThick_Dop_QuadInt_PG, numin, numax, args=(t[i], Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1,  full_output=fo  )[0] )	
	 return np.array(res) 



#########
###THICK
#########
def FThRnu_Thick_Dop_QuadInt(thet, r, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Thick_Dop, 0.,2.*ma.pi, args=(thet, r, nu, t, Dist, Aargs, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def FRnu_Thick_Dop_QuadInt(r, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThRnu_Thick_Dop_QuadInt, 0., ma.pi, args=(r, nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def Fnu_Thick_Dop_QuadInt(nu, t, Dist, Aargs, RHStable, Ttable):
	#Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	n0 = Aargs[5]
	Rd = Aargs[6]
	p  = Aargs[7]
	aeff = Aargs[10]
	Rtau1 = Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) + 0.01*Rd
	return intg.quad(FRnu_Thick_Dop_QuadInt, Rd, Rtau1, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def F_Thick_Dop_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnu_Thick_Dop_QuadInt, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

	# res=[]
	# for i in range(len(t)):
	# 	res = intg.quad(FRnu_Thick_Dop_QuadInt, numin, numax, args=(nu, t[i], Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	# return np.array(res)


def FThRnu_Thick_Dop_QuadInt_PG(thet, r, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(Fnuint_Thick_Dop_PG, 0.,2.*ma.pi, args=(thet, r, nu, t, Dist, Aargs, RHStable, Ttable) , epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def FRnu_Thick_Dop_QuadInt_PG(r, nu, t, Dist, Aargs, RHStable, Ttable):
	return intg.quad(FThRnu_Thick_Dop_QuadInt_PG, 0., ma.pi, args=(r, nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def Fnu_Thick_Dop_QuadInt_PG(nu, t, Dist, Aargs, RHStable, Ttable):
	#Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
	n0 = Aargs[5]
	Rd = Aargs[6]
	p  = Aargs[7]
	aeff = Aargs[10]
	Rtau1 = Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) + 0.01*Rd
	return intg.quad(FRnu_Thick_Dop_QuadInt_PG, Rd, Rtau1, args=(nu, t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]

def F_Thick_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable):
	#return intg.quad(Fnu_Thick_Dop_QuadInt_PG, numin, numax, args=(t, Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	#
	res=[]
	for i in range(len(t)):
		res.append( intg.quad(Fnu_Thick_Dop_QuadInt_PG, numin, numax, args=(t[i], Dist, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0] )
	return np.array(res)

































































####################################################
###  OLD INTEGRATION   
####################################################


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
	n0 = n0 * (pp-1.)/(ma.pi * aeff*aeff * Rde)   # must be greater than one
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
	
	n0 = Aargs[5]
	Rd = Aargs[6]
	p  = Aargs[7]
	aeff = Aargs[10]
	
	Rout = Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) + 0.01*Rd

	return intg.quad(Fnudphidth_Thick, Rd, Rout, args=(nu, t, Dist, Rout, Aargs, RHStable, Ttable), epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1 )[0]

def Fnu_Thick_mult(nu, t, Dist, Rout, Aargs, RHStable, Ttable):#,tauGrid):
	n0 = Aargs[5]
	Rd = Aargs[6]
	p  = Aargs[7]
	aeff = Aargs[10]
	
	Rout = Rd*(  1. - (p - 1.)/(n0*ma.pi*aeff*aeff*Rd)  )**(1./(1. - p)) + 0.01*Rd
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

