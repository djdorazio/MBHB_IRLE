import cPickle as pickle


import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt


matplotlib.rcParams.update({'font.size': 20})


from FluxFuncs_IRLE import *
from ErrFuncs_IRLE import *


import scipy.integrate as intgt



ISOvJ= False
DOPvJ= True





Nnumn = 0.0
Nnumx = 3.0

numn = Nnumn*numicron
numx = Nnumx*numicron


##ISO VARS:
Amp = 0.22105  ##For alpha=0.0, beta=0.0677
t0 = ma.pi/2. #0.0




#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
Omb = 1./(1*yr2sec)
L0 = 6.78*10**46
MPGmx = 10**9.4*Msun
Ryr = c*yr2sec
RdPG = 3.0*pc2cm#ma.sqrt(1.)*2.8 *pc2cm
OmPG = 2.*ma.pi/(1474*3600*24) #Omb*2.*ma.pi/4.1
Ombn = OmPG
#alphnu = 1.1

Rorb = c*2.*ma.pi/Omb
Ompc = 2.*ma.pi*c/pc2cm/2.


## TEST VALUES
### DUST stuff
## for Qv
nne = 0.0#1.8
nu0 = numicron#/0.37

Rde = RdPG
Rrout = 1.0*Rde
pp = 2.0
thetTst = 0.0#1.*ma.pi/2.


aeff = (c/nu0)/(2.*ma.pi)
#aeff = 0.16*10**(-4) #(0.1 micrometer is an average ISM dust grain size - choose 0.16 to make nu0~1um)
#md = 10**(-14)
n10 = 1.0/(ma.pi*Rde*aeff*aeff) * (pp-1.)
nfac = 10.0
n0 = nfac*n10 ##6.*10**5*Msun/md * 1./(4./3.*ma.pi*(Rrout**3 - Rde**3))


Lav = L0
betst = 0.06776  ## gives 0.14 mag amplitdue at I=0 
Inc = 0.0#ma.pi/2#ma.acos(0.07/betst)#0.*np.pi/4.
#Ombn = 2.*ma.pi/(Rde/c)  ##(2pi/P)
alph = 4.0
Dst = 1.4*10**9*pc2cm




### WISE BAND + Observational STUFF
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3


W1_mid = (W1mx + W1mn)/2.
W2_mid = (W2mx + W2mn)/2.

W3_mid = numicron/12.
W4_mid = numicron/22.



nuVbnd = c/(5.45*10**(-5))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-21)*8.8560*10**(13)#(W1mn + W1mx)/2
FW2Rel = 1.71787*10**(-21)*6.4451*10**(13)#(W2mn + W2mx)/2


# if (not (Qv_k or Qv_nu0)):
# 	##TABULATE T's and RHSs
# 	print "Creating look up tables"
# 	NT = 10000
# 	RHS_table = np.zeros(NT)
# 	T_table = np.linspace(1., 2000., NT)
# 	for i in range(NT):
# 		RHS_table[i] = T_RHS(T_table[i], nu0, nne)







Tmin = 100.
Tsub = 1800.
NT   = 1800



Rsub = 0.5 * ma.sqrt(Lav/(10.**(46))) * (1800./1800.)**(2.8)  * pc2cm

frc_PG = Rsub/c / (2.*ma.pi/OmPG) ## sublimation radius for L = Lav (10^8.7 is Mass for which epsilon = 1)

#frc_PG = Rde /ma.sqrt(0.1)  * ma.sqrt( 10**(8.7)/(10**(9.0)) )/c / (2.*ma.pi/OmPG) ## sublimation radisu for L = Lav (10^8.7 is Mass for which epsilon = 1)
#frc_PGb = frc_PG /ma.sqrt(0.1)

frc_mx = 3.0
frc_t = np.linspace(0.01, frc_mx, 40.0)
frc_a = np.linspace(0.01, frc_mx, 1000.0)

WeinCst = 2.821439
TW1_Wein = 1./ WeinCst * h * W1_mid/kb
TW2_Wein = 1./ WeinCst * h * W2_mid/kb


RW1_Wein = 0.5 * ma.sqrt(Lav/(10.**(46))) * (1800./TW1_Wein)**(2.8)  * pc2cm
RW2_Wein = 0.5 * ma.sqrt(Lav/(10.**(46))) * (1800./TW2_Wein)**(2.8)  * pc2cm

#RW1_Wein = ma.sqrt(Lav /( 16. * ma.pi * sigSB * TW1_Wein**4)  )
#RW2_Wein = ma.sqrt(Lav /( 16. * ma.pi * sigSB * TW2_Wein**4)  )

#RW1_Wein = ma.sqrt(4. * Lav /( 9. * ma.pi * sigSB * TW1_Wein**4)  )
#RW2_Wein = ma.sqrt(4. * Lav /( 9. * ma.pi * sigSB * TW2_Wein**4)  )

tdoP_W1_Wein = RW1_Wein/c/(2*ma.pi/OmPG) 
tdoP_W2_Wein = RW2_Wein/c/(2*ma.pi/OmPG) 



## First Compute amplitdue
def ISO_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc

	#Ombn = OmPG
	#Rd = 2. *ma.pi * c/OmPG * frc
	t0 = 0.0
	arg2 = [Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	LIR_mx = F_Sphere_Iso_QuadInt(numn, numx, 0.25*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Sphere_Iso_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff
	
	## INTEGRATE WITH SERIES EXPANSION ()gam = 0 for now
	#LIR_mx = F_Sphere_Iso_QuadInt_Series( 0.25*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table, numn, numx)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	#LIR_mn = F_Sphere_Iso_QuadInt_Series( 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table, numn, numx)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff
	

	#LIR_bol = F_Sphere_Iso_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	#return np.abs(0.5*(LIR_mx - LIR_mn)/(Lav*Amp))
	
	return 0.5*(LIR_mx - LIR_mn)/(Lav*Amp) *   Lav/(0.5*(LIR_mx + LIR_mn))
	#return np.log10(LIR_mx/LIR_mn)/np.log10( (1.+Amp)/(1.-Amp) )


## First Compute amplitdue
def ISO_MAX(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc
	t0 = 0.0
	arg2 = [Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	LIR_mx = F_Sphere_Iso_QuadInt(numn, numx, 0.25*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = 0.0#F_Sphere_Iso_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff
	

	#LIR_bol = F_Sphere_Iso_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	#return np.abs(0.5*(LIR_mx - LIR_mn)/(Lav*Amp))
	
	return 0.5*(LIR_mx - LIR_mn)/(Lav*Amp) #*   Lav/(0.5*(LIR_mx + LIR_mn))
	#return np.log10(LIR_mx/LIR_mn)/np.log10( (1.+Amp)/(1.-Amp) )




# For optically thick case
def ISO_OptThick_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc
	t0 = 0.0
	arg2 = [Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	LIR_mx = F_ShTorOptThick_Iso_QuadInt(numn, numx, 0.25*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_ShTorOptThick_Iso_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	        

	#LIR_bol = F_Sphere_Iso_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	#return np.abs(0.5*(LIR_mx - LIR_mn)/(Lav*Amp))
	return 0.5*(LIR_mx - LIR_mn)/(Lav*Amp) *   Lav/(0.5*(LIR_mx + LIR_mn))


def DOP_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc

	#Ombn = OmPG
	#Rd = 2. *ma.pi * c/OmPG * frc
	
	#nne = 0.0
	#thetTst = 
	arg2 = [Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	
	LUV_mx = Fsrc_Dop(0.5*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LUV_mn = Fsrc_Dop(0.75*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
	LIR_mx = F_Sphere_Dop_QuadInt(numn, numx, (0.25 +0.25)*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Sphere_Dop_QuadInt(numn, numx, (0.5+0.25)*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	

	#LIR_bol = F_Sphere_Dop_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 


	#return np.abs((LIR_mx - LIR_mn)/(LUV_mx - LUV_mn))
	return (LIR_mx - LIR_mn)/(LUV_mx - LUV_mn) * Lav/(0.5*(LIR_mx + LIR_mn))
	#return np.log10(LIR_mx/LIR_mn)/np.log10(LUV_mx/LUV_mn)

def DOP_OptThick_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc

	arg2 = [Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	
	LUV_mx = Fsrc_Dop(0.5*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LUV_mn = Fsrc_Dop(0.75*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
	LIR_mx = F_ShTorOptThick_Dop_QuadInt(numn, numx, 0.5*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_ShTorOptThick_Dop_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	

	#LIR_bol = F_Sphere_Dop_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 

	return (LIR_mx - LIR_mn)/(LUV_mx - LUV_mn) * Lav/(0.5*(LIR_mx + LIR_mn))




def DOP_AIR_o_AUV_NRM(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc
	#nne = 0.0
	#thetTst = 
	arg2 = [Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	
	LUV_mx = Fsrc_Dop(0.5*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LUV_mn = Fsrc_Dop(0.75*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
	LIR_mx = F_Sphere_Dop_QuadInt(numn, numx, 0.5*2*ma.pi/Ombn, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Sphere_Dop_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
	#return np.abs((LIR_mx - LIR_mn)/(LUV_mx - LUV_mn))
	AIR = np.exp( (LIR_mx - LIR_mn)*0.5/Lav )
	AUV = np.exp( (LUV_mx - LUV_mn)*0.5/Lav )

	AIR = (LIR_mx - LIR_mn)*0.5/Lav 
	AUV = (LUV_mx - LUV_mn)*0.5/Lav 

	return AIR/AUV




































if (ISOvJ):
	thetTst = ma.pi/2.05

	J1 = 0.0
	J2 = ma.pi/8.
	J3 = ma.pi/4.
	J4 = 3.*ma.pi/8.
	J5 = ma.pi/2.0
	
	

	nne = 0.0 #no absorption efficiency
	##TABULATE T's and RHSs
	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)


	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, J1, aeff, nu0, nne]
	arg2_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, J2, aeff, nu0, nne]
	arg3_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, J3, aeff, nu0, nne]
	arg4_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, J4, aeff, nu0, nne]
	arg5_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, J5, aeff, nu0, nne]


	from scipy import special as spc
	J0_anl   =  spc.j0(2.*ma.pi*frc_a)
	Jpo2_anl =  1./(2.*ma.pi*frc_a* ma.cos(thetTst)) * np.sin(2.*ma.pi*frc_a * ma.cos(thetTst)) 
	


	
	ISO1_AIR_over_AUV = np.zeros(len(frc_t))
	ISO2_AIR_over_AUV = np.zeros(len(frc_t))
	ISO3_AIR_over_AUV = np.zeros(len(frc_t))
	ISO4_AIR_over_AUV = np.zeros(len(frc_t))
	ISO5_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		ISO1_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
		ISO2_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_ISO, RHS_table, T_table)
		ISO3_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_ISO, RHS_table, T_table)
		ISO4_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg4_ISO, RHS_table, T_table)
		ISO5_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg5_ISO, RHS_table, T_table)




	plt.figure()

	plt.title(r'Isotropic, $\theta_T = \pi/2.05$')
	Jpo2_ANL  = plt.plot(frc_a, Jpo2_anl, color='black',  linewidth=3)
	J0_ANL    = plt.plot(frc_a, J0_anl, color='green',  linewidth=3)

	

	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='green', s=40)
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='blue', s=40)
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='red', s=40)
	ISO4 = plt.scatter(frc_t, ISO4_AIR_over_AUV, marker='x', color='purple', s=40)
	ISO5 = plt.scatter(frc_t, ISO5_AIR_over_AUV, marker='x', color='black', s=40)
	plt.plot(frc_t, ISO1_AIR_over_AUV, linestyle='--', color='green')
	plt.plot(frc_t, ISO2_AIR_over_AUV, linestyle='--', color='blue')
	plt.plot(frc_t, ISO3_AIR_over_AUV, linestyle='--', color='red')
	plt.plot(frc_t, ISO4_AIR_over_AUV, linestyle='--', color='purple')
	plt.plot(frc_t, ISO5_AIR_over_AUV, linestyle='--', color='black')


	plt.axhline(y=0, color='black', linestyle=':')


	plt.legend( [ J0_ANL[0], ISO1, ISO2, ISO3, ISO4, Jpo2_ANL[0], ISO5], (r'$J = 0$', '',  r'$J = \pi/8$',    r'$J = \pi/4$',  r'$J = 3\pi/8$', r'$J = \pi/2$', ''), loc='upper right', fontsize=18)


	plt.axhline(y=0, color='black', linestyle=':')



	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")



	plt.xlim(0.0,frc_mx)
	plt.ylim(-0.7,1.0)

	plt.tight_layout()
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/AIRoAUV_Jdep_thT%g_numin%g_numx%g_reclim1.png" %(thetTst, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)








if (DOPvJ):
	thetTst = ma.pi/2.05

	J1 = 0.0
	J2 = ma.pi/8.
	J3 = ma.pi/4.
	J4 = 3.*ma.pi/8.
	J5 = ma.pi/2.0

	#ni = 1./4.
	Inc = ma.pi/2.
	
	

	nne = 0.0 #no absorption efficiency
	##TABULATE T's and RHSs
	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)


	arg1_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, J1, aeff, nu0, nne]
	arg2_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, J2, aeff, nu0, nne]
	arg3_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, J3, aeff, nu0, nne]
	arg4_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, J4, aeff, nu0, nne]
	arg5_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, J5, aeff, nu0, nne]




	Dop1_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst))/np.cos(thetTst) )

	from scipy import special as spc
	DopRing_anal = spc.j1(2.*ma.pi*frc_a)

	DOP1_AIR_over_AUV = np.zeros(len(frc_t))
	DOP2_AIR_over_AUV = np.zeros(len(frc_t))
	DOP3_AIR_over_AUV = np.zeros(len(frc_t))
	DOP4_AIR_over_AUV = np.zeros(len(frc_t))
	DOP5_AIR_over_AUV = np.zeros(len(frc_t))


	for i in range(len(frc_t)):
		DOP1_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		DOP2_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_DOP, RHS_table, T_table)
		DOP3_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_DOP, RHS_table, T_table)
		DOP4_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg4_DOP, RHS_table, T_table)
		DOP5_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg5_DOP, RHS_table, T_table)



	plt.figure()


	plt.title(r'Doppler, $\theta_T = \pi/2.05$, $I = \pi/2$') #$I = %g$' %Inc)


	Anl1 = plt.plot(frc_a, Dop1_anal, color='black',  linewidth=3)
	Anl2 = plt.plot(frc_a, DopRing_anal, color='green',  linewidth=3)

	DOP1 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='green', s=40)
	DOP2 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='blue', s=40)
	DOP3 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='red', s=40)
	DOP4 = plt.scatter(frc_t, DOP4_AIR_over_AUV, marker='*', color='purple', s=40)
	DOP5 = plt.scatter(frc_t, DOP5_AIR_over_AUV, marker='*', color='black', s=40)
	plt.plot(frc_t, DOP1_AIR_over_AUV, color='green', linestyle="--")
	plt.plot(frc_t, DOP2_AIR_over_AUV, color='blue', linestyle="--")
	plt.plot(frc_t, DOP3_AIR_over_AUV, color='red', linestyle="--")
	plt.plot(frc_t, DOP4_AIR_over_AUV, color='purple', linestyle="--")
	plt.plot(frc_t, DOP5_AIR_over_AUV, color='black', linestyle="--")


	plt.axhline(y=0, color='black', linestyle=':')


	plt.legend( [ Anl2[0],  DOP1, DOP2, DOP3, DOP4, Anl1[0], DOP5 ], (r'$J = 0$', '',  r'$J = \pi/8$',    r'$J = \pi/4$',  r'$J = 3\pi/8$', r'$J = \pi/2$', ''), loc='upper right', fontsize=18)


	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")



	plt.xlim(0.0,frc_mx)
	plt.ylim(-0.5,1.0)

	plt.tight_layout()
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/AIRoAUVDop_Jdep_thT%g_Inc%gnumin%g_numx%g_reclim2.png" %(thetTst, Inc, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)














