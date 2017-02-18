import cPickle as pickle


import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt


matplotlib.rcParams.update({'font.size': 16})


from FluxFuncs_IRLE import *
from ErrFuncs_IRLE import *


import scipy.integrate as intgt

from scipy import special as spc


#push test
###OPTIONS###OPTIONS
IR_Lum = False
AddW4 = False #W1, W2, W3, W4 if True, otherwise see AddW3
AddW3 = True #W1, W2, W3 if True, #W1, W2 if False
fitQv = False
fitFcov = True
func_min = False
MCMC = True
NThread = 8

##FIG 3
ISOvDop = False

ISOvDop_MAX = False

##Fig 4
Dop_alphs = False
Dop_alphs_Rings = False

ISOvDop_varyI = False

#Figs 8 and 10
JJt = ma.pi/2.
ISOvthT = False
DOPvthT = False

PG1302_ISO = False
PG1302_ISO_Jth_Conts = False
PG1302_Dop = False
PG1302_Dop_Jth_Conts = True


W1W2vals = True
numRing = True


##Fig 5
Qv_nu0 = False
## INTEGRATION LIMTS FOR ALL nu
if (Qv_nu0==True and ISOvDop == False and Dop_alphs == False and ISOvthT == False and DOPvthT == False and numRing == False):
	print "limiting nu"
	Nnumn = 1./2.8 #0.0#0.00001
	Nnumx = 1./4.0 #5.0
else:
	Nnumn = 0.0
	Nnumx = 3.0

numn = Nnumn*numicron
numx = Nnumx*numicron

Qv_k = False


ThkvThn_ISO = False
ThkvThn_Dop = False
fit_Dop = False





##ISO VARS:
Amp = 0.15# For alph = 1, bet = 0.07 0.22105  ##For alpha=0.0, beta=0.0677
t0 = ma.pi/2. #0.0




#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
Omb = 1./(1*yr2sec)
L0 = 6.78*10**46
MPGmx = 10**9.4*Msun
Ryr = c*yr2sec
RdPG = 1.0*pc2cm#1.0*pc2cm#ma.sqrt(1.)*2.8 *pc2cm
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
betst = 0.07  ## gives 0.14 mag amplitdue at I=0, alpha=1.1
Inc = 0.0#ma.pi/2#ma.acos(0.07/betst)#0.*np.pi/4.
#Ombn = 2.*ma.pi/(Rde/c)  ##(2pi/P)
alph = 1.0  ## for getting optical LC
Dopalph = -1.0 ## for illuminate dust
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







Tmin = 10.
Tsub = 1800.
NT   = 1890



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


##IR LUM
def IRLum(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc #don't worry beta isn't set by Ombn
	t0 = 0.0
	arg2 = [Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	LIR_mx = F_Sphere_Iso_QuadInt(numn, numx, 0.25*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Sphere_Iso_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff
	

	#LIR_bol = F_Sphere_Iso_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	#return np.abs(0.5*(LIR_mx - LIR_mn)/(Lav*Amp))
	
	return 0.5*(LIR_mx + LIR_mn)

## First Compute amplitdue
def ISO_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc #don't worry beta isn't set by Ombn

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




def ISO_Ring_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1

	Ombn = 2. *ma.pi * c/Rd * frc #don't worry beta isn't set by Ombn

	t0 = 0.0
	arg2 = [Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	LIR_mx = F_Ring_Iso_QuadInt( 0.25*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Ring_Iso_QuadInt( 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff
	
	
	return 0.5*(LIR_mx - LIR_mn)/(Lav*Amp) *   Lav/(0.5*(LIR_mx + LIR_mn))




## First Compute amplitdue
def ISO_MAX(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc #don't worry beta isn't set by Ombn
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
	Ombn = 2. *ma.pi * c/Rd * frc #don't worry beta isn't set by Ombn
	t0 = 0.0
	arg2 = [Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	LIR_mx = F_ShTorOptThick_Iso_QuadInt(numn, numx, 0.25*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_ShTorOptThick_Iso_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	        

	#LIR_bol = F_Sphere_Iso_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	#return np.abs(0.5*(LIR_mx - LIR_mn)/(Lav*Amp))
	return 0.5*(LIR_mx - LIR_mn)/(Lav*Amp) *   Lav/(0.5*(LIR_mx + LIR_mn))


def DOP_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc  #don't worry beta isn't set by Ombn


	#Ombn = OmPG
	#Rd = 2. *ma.pi * c/OmPG * frc
	
	#nne = 0.0
	#thetTst = 
	arg2 = [Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	
	#LUV_mx = Fsrc_Dop(0.5*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	#LUV_mn = Fsrc_Dop(0.75*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
	LUV_mx = Fsrc_Dop(0.25*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LUV_mn = Fsrc_Dop(0.5*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	

	LIR_mx = F_Sphere_Dop_QuadInt(numn, numx, (0.25-0.25)*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Sphere_Dop_QuadInt(numn, numx, (0.5-0.25)*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	

	#LIR_bol = F_Sphere_Dop_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 


	#return np.abs((LIR_mx - LIR_mn)/(LUV_mx - LUV_mn))
	return (LIR_mx - LIR_mn)/(LUV_mx - LUV_mn) * Lav/(0.5*(LIR_mx + LIR_mn))
	#return np.log10(LIR_mx/LIR_mn)/np.log10(LUV_mx/LUV_mn)




def DOP_Ring_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc

	arg2 = [Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	
	LUV_mx = Fsrc_Dop(0.5*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LUV_mn = Fsrc_Dop(0.75*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
	LIR_mx = F_Ring_Dop_QuadInt( (0.25 +0.25)*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Ring_Dop_QuadInt( (0.5+0.25)*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 

	return (LIR_mx - LIR_mn)/(LUV_mx - LUV_mn) * Lav/(0.5*(LIR_mx + LIR_mn))







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


def DOP_MAX(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc
	#nne = 0.0
	#thetTst = 
	arg2 = [Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	
	LUV_mx = Fsrc_Dop(0.5*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LUV_mn = Fsrc_Dop(0.75*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
	LIR_mx = F_Sphere_Dop_QuadInt(numn, numx, 0.5*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = 0.0#F_Sphere_Dop_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	

	#LIR_bol = F_Sphere_Dop_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 


	#return np.abs((LIR_mx - LIR_mn)/(LUV_mx - LUV_mn))
	return (LIR_mx - LIR_mn)/(LUV_mx - LUV_mn) #* Lav/(0.5*(LIR_mx + LIR_mn))
	#return np.log10(LIR_mx/LIR_mn)/np.log10(LUV_mx/LUV_mn)


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



def Dop_Anyl_LC(t,frc, Dst, arg1):
	Lav, bet, Inc, Om, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	gam = 1./np.sqrt(1.-bet*bet)
	td = Rd/c
	return 	(gam*(64*Om**3*td**3 + 128*bet**2*Om**3*td**3 - 36*bet**2*Om*td*np.cos(2*Om*t) - 18*bet**2*Om*td*np.cos(2*(Inc - Om*t)) - 18*bet**2*Om*td*np.cos(2*(Inc + Om*t)) - 36*bet**2*Om*td*np.cos(2*Om*(t - 2*td)) + 256*bet*Om**2*td**2*np.cos(Inc)*np.cos(Om*(t - td))*np.cos(Om*td) - 18*bet**2*Om*td*np.cos(2*(Inc + Om*t - 2*Om*td)) - 18*bet**2*Om*td*np.cos(2*(Inc - Om*t + 2*Om*td)) + 18*bet**2*np.sin(2*Om*t) - 18*bet**2*Om**2*td**2*np.sin(2*Om*t) - 9*bet**2*np.sin(2*(Inc - Om*t)) + 9*bet**2*np.sin(2*(Inc + Om*t)) - 9*bet**2*Om**2*td**2*np.sin(2*(Inc + Om*t)) - 18*bet**2*np.sin(2*Om*(t - 2*td)) + 18*bet**2*Om**2*td**2*np.sin(2*Om*(t - 2*td)) - 256*bet*Om*td*np.cos(Inc)*np.cos(Om*(t - td))*np.sin(Om*td) - 12*bet**2*Om**2*td**2*np.cos(2*Om*(t - td))*np.sin(2*Om*td) - 6*bet**2*Om**2*td**2*np.cos(2*(Inc + Om*t - Om*td))*np.sin(2*Om*td) - 24*bet**2*Om**2*td**2*np.cos(2*(Inc - Om*t + Om*td))*np.sin(2*Om*td) - 9*bet**2*np.sin(2*(Inc + Om*t - 2*Om*td)) + 9*bet**2*Om**2*td**2*np.sin(2*(Inc + Om*t - 2*Om*td)) + 9*bet**2*np.sin(2*(Inc - Om*t + 2*Om*td))))/(64.*Om**3*td**3)





# def DOP_Anyl_AIRoA(frc, numn, numx, Dst, arg1, RHS_table, T_table):
# 	Lav, bet, Inc, Om, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1

# 	Om = 2. *ma.pi * c/Rd * frc

# 	arg2 = [Lav, bet, Inc, Om, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]


	
# 	alph = -1.0


# 	LUV_mx = Fsrc_Dop(0.5*2*ma.pi/Om, Dst, ma.pi/2., 0.0, Lav, bet, Inc, Om, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
# 	LUV_mn = Fsrc_Dop(0.75*2*ma.pi/Om, Dst, ma.pi/2., 0.0, Lav, bet, Inc, Om, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
# 	LIR_mx = Dop_Anyl_LC((0.5)*2*ma.pi/Om + Rd/c, frc, Dst, arg2)
# 	LIR_mn = Dop_Anyl_LC((0.75)*2*ma.pi/Om + Rd/c, frc, Dst, arg2)


# 	return (LIR_mx - LIR_mn)/(LUV_mx - LUV_mn) * Lav/(0.5*(LIR_mx + LIR_mn))



def DOP_Anyl_AIRoA(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, bet, Inc, Om, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1

	td = Rd/c
	Om = 2. *ma.pi /td * frc

	return -(6*(-1 + bet*np.cos(Inc))**4*(2*Om*td*(np.cos(thetTst)*(8*Om*td*np.cos(Om*td*np.cos(thetTst)) + 15*bet*np.cos(Inc)*np.cos(2*Om*td*np.cos(thetTst))) - 8*np.sin(Om*td*np.cos(thetTst))) + 5*bet*np.cos(Inc)*(-3 + Om**2*td**2 + 3*Om**2*td**2*np.cos(2*thetTst))*np.sin(2*Om*td*np.cos(thetTst))))/((-1 + bet**2)**2*Om*td*(-2 + bet*np.cos(Inc))*(2 + bet*np.cos(Inc)*(-2 + bet*np.cos(Inc)))*(Om*td*np.cos(thetTst)*(-((24 + 85*bet**2)*Om*td) + bet*(48*np.cos(Inc)*np.cos(Om*td*np.cos(thetTst)) + 5*bet*Om*td*(np.cos(2*thetTst) + 6*np.cos(2*Inc)*np.sin(thetTst)**2))) - 48*bet*np.cos(Inc)*np.sin(Om*td*np.cos(thetTst))))




def DOP_EORing2_Anyl_AIRoA(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, bet, Inc, Om, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1

	td = Rd/c
	Om = 2. *ma.pi /td * frc

	return -( (4*spc.j1(Om*td)*(1 + bet*np.cos(Inc))**4)/((-1. + bet**2)**2*(-1 + 2*bet*spc.j1(Om*td)*np.cos(Inc))*(4 + 6*bet*np.cos(Inc) + 4*bet**2*np.cos(Inc)**2 + bet**3*np.cos(Inc)**3)) )





def DOP_EORing3_Anyl_AIRoA(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, bet, Inc, Om, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1

	td = Rd/c
	Om = 2. *ma.pi /td * frc

	return -( (-2*(1 + bet*np.cos(Inc))**4*(4*Om*td*np.cos(Inc)*(2*spc.j1(Om*td) + 5*bet*spc.j0(2*Om*td)*np.cos(Inc)) - 5*bet*spc.j1(2*Om*td)*(3 + np.cos(2*Inc)))*  1./np.cos(Inc)   )/((-1 + bet**2)**2*Om*td*(4 + 6*bet*np.cos(Inc) + 4*bet**2*np.cos(Inc)**2 + bet**3*np.cos(Inc)**3)*(4 + 15*bet**2 - 8*bet*spc.j1(Om*td)*np.cos(Inc) + 5*bet**2*np.cos(2*Inc))) )



























#### BB FIT BELOW
if (IR_Lum):
	# nne = 0.0 #no absorption efficiency
	# nu0 = numicron
	# ##TABULATE T's and RHSs
	# print "Creating look up tables"
	# Tsub = 1800.
	# NT   = 1800
	# RHS_table = np.zeros(NT)
	# T_table = np.linspace(0.01, Tsub, NT)
	# for i in range(NT):
	# 	RHS_table[i] = T_RHS(T_table[i], nu0, nne)
	# arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]

	# LIR = IRLum(1., numn, numx, Dst, arg1_ISO, RHS_table, T_table)
	# LW1 = FW1Rel * 10.**(-(11.3)/2.5) *(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	# LW2 = FW2Rel * 10.**(-(10.25)/2.5) *(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 

	# print "LIR(sphere) = LUV = %g" %LIR
	# print "LW1 = %g" %LW1
	# print "LW2 = %g" %LW2


	#read in fluxes and errors, take average
	FnuW1_t = np.genfromtxt("../dat/W1.txt",usecols=0) * 10.**(-26)
	FnuW2_t = np.genfromtxt("../dat/W2.txt",usecols=0) * 10.**(-26)
	ErrW1_t = np.genfromtxt("../dat/W1.txt",usecols=1) * 10.**(-26)
	ErrW2_t = np.genfromtxt("../dat/W2.txt",usecols=1) * 10.**(-26)


	Fw1 = np.mean(FnuW1_t) 
	Fw2 = np.mean(FnuW2_t)
	Fw3 = 36.0745 * 10.**(-26)
	Fw4 = 101.8591 * 10.**(-26)


	ErrW1 = np.mean(ErrW1_t)
	ErrW2 = np.mean(ErrW2_t)
	ErrW3 = 0.5647 * 10.**(-26)
	ErrW4 = 2.6262 * 10.**(-26)


	## BEST FIT Black Body
	from scipy import optimize
	from scipy.optimize import fmin

	from scipy import special as spc
	from emcee_Funcs import *

	if (AddW4):
		# All
		Shell_File = "AddW4"
		nus = (W1_mid, W2_mid, W3_mid, W4_mid)
		nus = np.array(nus)
		Flxs = (Fw1, Fw2, Fw3, Fw4)
		Flxs = np.array(Flxs)
		Errs = (ErrW1, ErrW2, ErrW3, ErrW4)
		Errs = np.array(Errs)
	elif (AddW3):
		# NO W4
		Shell_File = "AddW3"
		nus = (W1_mid, W2_mid, W3_mid)
		nus = np.array(nus)
		Flxs = (Fw1, Fw2, Fw3)
		Flxs = np.array(Flxs)
		Errs = (ErrW1, ErrW2, ErrW3)
		Errs = np.array(Errs)
	else:
		#No W# and no W4
		Shell_File = "W1W2"
		nus = (W1_mid, W2_mid)
		nus = np.array(nus)
		Flxs = (Fw1, Fw2)
		Flxs = np.array(Flxs)
		Errs = (ErrW1, ErrW2)
		Errs = np.array(Errs)	

	

	#ALL
	# nus = (W1_mid, W2_mid, W3_mid, W4_mid)
	# nus = np.array(nus)
	# Flxs = (Fw1, Fw2, Fw3, Fw4)
	# Flxs = np.array(Flxs)
	# Errs = (ErrW1, ErrW2, ErrW3, ErrW4)
	# Errs = np.array(Errs)

	if (AddW3):
		T0 = 570.0#0.5*(TW2_Wein + TW1_Wein)
	else:
		T0 = 880.0
	if (func_min):
		if (fitQv):
			Topt  = sc.optimize.fmin(BB_Err2_Qv,    [T0, 1.*numicron/10**14, 1.0, 1.0], args=(nus, Flxs, Errs), full_output=1, disp=False)[0]
			Td = Topt[0]
			nu0 = Topt[1] * 10**14
			gam = Topt[2]
			sqtfR = Topt[3] * pc2cm
			#cf1 = Topt[3]
		else:
			Topt  = sc.optimize.fmin(BB_Err2,    [T0, 1.0], args=(nus, Flxs, Errs), full_output=1, disp=False)[0]
			Td = Topt[0]
			#nu0 = numicron/0.37 defined above
			gam = nne
			sqtfR = Topt[1] * pc2cm
	else: 
		print "Not doing fmin"


	if (MCMC):
		import emcee

		if (fitQv):
			if (fitFcov):
				param_names = [r"$T_d$", r"$\nu_0$", r"$k$", r"$\cos{\theta_T}$", "BC/10"]
				BB_p0 = [T0, 0.1*numicron/10**14, 0.1, 0.125, 1.0]
			else:
				param_names = [r"$T_d$", r"$\nu_0$", r"$k$", r"$\sqrt{\cos{\theta_T}} R_d$"]
				if (func_min):
					BB_p0 = [Topt[0], Topt[1], Topt[2], Topt[3]]
				else:
					BB_p0 = [T0, 1.*numicron/10**14, 1.0, 0.125]
		else:	
			if (fitFcov):
				param_names = [r"$T_d$", r"$\cos{\theta_T}$", "BC/10"]
				BB_p0 = [T0, 0.1, 1.0]	
			else:
				param_names = [r"$T_d$", r"$\sqrt{\cos{\theta_T}} R_d$"]
				if (func_min):
					BB_p0 = [Topt[0], Topt[1]]
				else:
					BB_p0 = [T0, 1.0]		

		ndim = len(BB_p0)
		nwalkers = ndim*16#32

		if (fitQv):
			if (fitFcov):
				BB_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_BBposterior_Qv_Fcov, threads=NThread, args=(nus, Flxs, Errs) )
			else:
				BB_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_BBposterior_Qv, threads=NThread, args=(nus, Flxs, Errs) )
		else:
			if (fitFcov):
				BB_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_BBposterior_Fcov, threads=NThread, args=(nus, Flxs, Errs) )
			else:
				BB_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_BBposterior, threads=NThread, args=(nus, Flxs, Errs) )



		BB_p0 = np.array(BB_p0)
		BB_walker_p0 = np.random.normal(BB_p0, np.abs(BB_p0)*1E-4, size=(nwalkers, ndim))

		clen = 4096
		BB_pos,_,_ = BB_sampler.run_mcmc(BB_walker_p0 , clen)


		print "SAVING THE PICKLE mmmmm"
		with open("../emcee_data/Pickles/PG1302_BBfit_%iwalkers.pickle"%clen+Shell_File, "w") as f1:
			pickle.dump((BB_sampler.chain, BB_sampler.lnprobability), f1)

		### OPEN OUTPUT DATA
		with open("../emcee_data/Pickles/PG1302_BBfit_%iwalkers.pickle"%clen+Shell_File) as f1:
			BB_chain,BB_lnprobs = pickle.load(f1)


		BB_flatchain = np.vstack(BB_chain[:,clen/4:])
		BB_flatlnprobs = np.vstack(BB_lnprobs[:,clen/4:])

					
		BB_p_opt  = BB_flatchain[BB_flatlnprobs.argmax()]

		print "ANALYSING MCMC (!)..."
		with open("../emcee_data/Pickles/PG1302_BBfit_%iwalkers.pickle"%clen+Shell_File) as f1:
			BB_chain,BB_lnprobs = pickle.load(f1)


		##PLOT dem WALKERS
		for k in range(BB_chain.shape[2]):
			plt.figure()
			for i in range(BB_chain.shape[0]):
				plt.plot(BB_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
				plt.ylabel(param_names[k])
				plt.xlabel('steps')
				#print param_names[k]
				plt.tight_layout()
			plt.savefig('../emcee_data/BBfit_PG1302_%s_%iwalkers'%(param_names[k],clen)+Shell_File+'.png')
			plt.clf()


		###CORNER PLOT	
		BB_flatchain = np.vstack(BB_chain[:,clen/4:])
		BB_flatlnprobs = np.vstack(BB_lnprobs[:,clen/4:])

		#import triangle
		import corner as triangle
		BB_fig = triangle.corner(BB_flatchain, labels=param_names, quantiles=[0.15, 0.5, 0.85],show_titles=True, title_kwargs={"fontsize": 14},label_kwargs={"fontsize": 20})	
		BB_fig.savefig('../emcee_data/BBfit_PG1302_Corner_Plot_%iwalkers'%clen+Shell_File+'.png')


		## Do some stats on the walkers
		from scipy.stats import scoreatpercentile as scoretpercentile
		## max posterior + percentiles
		BB_MAP_vals = BB_flatchain[BB_flatlnprobs.argmax()]
		BB_perc = scoretpercentile(BB_flatchain, [15,85], axis=0)


		filename = "BBfit_resuts_%iwalkers"%clen+Shell_File+".txt" 
		print "Printing Results"
		target = open(filename, 'w')
		target.truncate()


		for i,name in enumerate(param_names):
			BB_diff_minus = BB_MAP_vals[i] - BB_perc[0,i]
			BB_diff_plus = BB_perc[1,i] - BB_MAP_vals[i]
			target.write("BBfit: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(BB_MAP_vals[i], BB_diff_plus, BB_diff_minus, name=name))
			target.write("\n")





		BB_mxprbs = zeros(nwalkers)
		
					

		for i in range(nwalkers):
			BB_mxprbs[i] = max(BB_lnprobs[i])
		

		chi2_pdf_BB = -max(BB_mxprbs)#/(len(nus) - len(param_names) - 1)
		
		target.write("\n")		
		target.write("BBfit reduced chi2 =  %04g" %chi2_pdf_BB)

		target.close()

		if (fitQv):
			Td  = BB_p_opt[0]
			nu0 = BB_p_opt[1] * 10**14
			gam = BB_p_opt[2] 
			
			if (fitFcov):
				Fcov = BB_p_opt[3]
				Lfac = BB_p_opt[4]
				Lav   = Lfac*Lav  #update L with Lfac then define R
				qIR   = (1./nu0)**(gam)
				RR    = ma.sqrt(  Lav / (4. * ma.pi * 8. * ma.pi  * qIR * h/c/c * (kb/h)**(4+gam) * spc.gamma(4+gam) * (spc.zetac(4+gam)+1.) * Td**(4+gam) ) )
				sqtfR = np.sqrt(Fcov) * RR

				### get avg deltas
				delT   = 0.5* ( (BB_MAP_vals[0] - BB_perc[0,0]) + (BB_perc[1,0] - BB_MAP_vals[0]) )
				delnu0 = 0.5* ( (BB_MAP_vals[1] - BB_perc[0,1]) + (BB_perc[1,1] - BB_MAP_vals[1]) ) * 10**14
				delk   = 0.5* ( (BB_MAP_vals[2] - BB_perc[0,2]) + (BB_perc[1,2] - BB_MAP_vals[2]) )
				delF    = 0.5* ( (BB_MAP_vals[3] - BB_perc[0,3]) + (BB_perc[1,3] - BB_MAP_vals[3]) ) 
				delLfac = 0.5* ( (BB_MAP_vals[4] - BB_perc[0,4]) + (BB_perc[1,4] - BB_MAP_vals[4]) ) 


				delR = ma.sqrt(  (1/2. *RR * delLfac/Lfac)**2  +  (1/2.* RR * delT/Td)**2  )
				delfR2 = ma.sqrt(  (RR*RR*delF)**2  +  (2. * Fcov * RR * delR)**2  )

				dfM = (BB_MAP_vals[3] - BB_perc[0,3])
				dfP = (BB_perc[1,3] - BB_MAP_vals[3])

				dLM = (BB_MAP_vals[3] - BB_perc[0,3])
				dLP = (BB_perc[1,3] - BB_MAP_vals[3])

			else:
				sqtfR = BB_p_opt[3] * pc2cm
				qIR   = (1./nu0)**(gam)
				RR = ma.sqrt(  Lav / (4. * ma.pi * 8. * ma.pi  * qIR * h/c/c * (kb/h)**(4+gam) * spc.gamma(4+gam) * (spc.zetac(4+gam)+1.) * Td**(4+gam) ) )
				Rg4 = ma.sqrt( Lav / (16.*ma.pi * sigSB*Td**4))


				### get avg deltas
				delT   = 0.5* ( (BB_MAP_vals[0] - BB_perc[0,0]) + (BB_perc[1,0] - BB_MAP_vals[0]) )
				delnu0 = 0.5* ( (BB_MAP_vals[1] - BB_perc[0,1]) + (BB_perc[1,1] - BB_MAP_vals[1]) ) * 10**14
				delk   = 0.5* ( (BB_MAP_vals[2] - BB_perc[0,2]) + (BB_perc[1,2] - BB_MAP_vals[2]) )
				delfR2 = 0.5* ( (BB_MAP_vals[3] - BB_perc[0,3]) + (BB_perc[1,3] - BB_MAP_vals[3]) ) * pc2cm
		else:
			if (fitFcov):
				Td = BB_p_opt[0]
				nu0 = 1.0
				gam = 0.0

				Fcov = BB_p_opt[1]
				Lfac = BB_p_opt[2]
				Lav   = Lfac*Lav  #update L with Lfac then define R
				qIR   = (1./nu0)**(gam)
				RR    = ma.sqrt(  Lav / (4. * ma.pi * 8. * ma.pi  * qIR * h/c/c * (kb/h)**(4+gam) * spc.gamma(4+gam) * (spc.zetac(4+gam)+1.) * Td**(4+gam) ) )
				sqtfR = np.sqrt(Fcov) * RR

				### get avg deltas
				## BB_perc[0,i] is lower value of ith parameter, [1,i] is higher value
				delT    = 0.5* ( (BB_MAP_vals[0] - BB_perc[0,0]) + (BB_perc[1,0] - BB_MAP_vals[0]) )
				delnu0  = 0.0
				delk    = 0.0
				delF    = 0.5* ( (BB_MAP_vals[1] - BB_perc[0,1]) + (BB_perc[1,1] - BB_MAP_vals[1]) ) 
				delLfac = 0.5* ( (BB_MAP_vals[2] - BB_perc[0,2]) + (BB_perc[1,2] - BB_MAP_vals[2]) ) 

				dfM = (BB_MAP_vals[1] - BB_perc[0,1])
				dfP = (BB_perc[1,1] - BB_MAP_vals[1])

				dLM = (BB_MAP_vals[2] - BB_perc[0,2])
				dLP = (BB_perc[1,2] - BB_MAP_vals[2])


				delR = ma.sqrt(  (1/2.* RR * delLfac/Lfac)**2  +  (1/2. *RR * delT/Td)**2  )
				delfR2 = ma.sqrt(  (RR*RR*delF)**2  +  (2. * Fcov * RR * delR)**2  )


			else:
				Td = BB_p_opt[0]
				#nu0 = numicron/0.37 defined above
				gam = nne
				sqtfR = BB_p_opt[1] * pc2cm

				qIR   = (1./nu0)**(gam)
				RR = ma.sqrt(  Lav / (4. * ma.pi * 8. * ma.pi  * qIR * h/c/c * (kb/h)**(4+gam) * spc.gamma(4+gam) * (spc.zetac(4+gam)+1.) * Td**(4+gam) ) )
				Rg4 = ma.sqrt( Lav / (16.*ma.pi * sigSB*Td**4))


				### get avg deltas
				delT   = 0.5* ( (BB_MAP_vals[0] - BB_perc[0,0]) + (BB_perc[1,0] - BB_MAP_vals[0]) )
				delnu0 = 0.0
				delk   = 0.0
				delfR2 = 0.5* ( (BB_MAP_vals[1] - BB_perc[0,1]) + (BB_perc[1,1] - BB_MAP_vals[1]) ) * pc2cm

	else:
		print "Not doing MCMC"


	

	







	#CFmin = 4.*ma.pi*Dst*Dst / Lav * (Fw2-ErrW1)/(ma.pi*Bv(W2_mid, Topt[0]))  * sigSB  * Topt[0]**4
	#CFmax = 4.*ma.pi*Dst*Dst / Lav * (Fw2+ErrW1)/(ma.pi*Bv(W2_mid, Topt[0]))  * sigSB  * Topt[0]**4

	# CF = 4.*ma.pi*(Rmatch*pc2cm)**2 * sigSB  * Topt[0]**4 / Lav
	# CFmin = 4.*ma.pi*(Rmatch*pc2cm)**2 * sigSB  * Topt[0]**4 / (Lav*15./10.)
	# CFmax = 4.*ma.pi*(Rmatch*pc2cm)**2 * sigSB  * Topt[0]**4 / (Lav*5./10.)

	
	# CF    = cf1*4.*ma.pi*RR**2 * (8. * ma.pi  * qIR * h/c/c * (kb/h)**(4.+gam) * spc.gamma(4.+gam) * spc.zetac(4.+gam) * Td**(4.+gam) ) / Lav
	# CFmin = cf1*4.*ma.pi*RR**2 * (8. * ma.pi  * qIR * h/c/c * (kb/h)**(4.+gam) * spc.gamma(4.+gam) * spc.zetac(4.+gam) * Td**(4.+gam) ) / (Lav*15./10.)
	# CFmax = cf1*4.*ma.pi*RR**2 * (8. * ma.pi  * qIR * h/c/c * (kb/h)**(4.+gam) * spc.gamma(4.+gam) * spc.zetac(4.+gam) * Td**(4.+gam) ) /(Lav*5./10.)



	#PROPAGATE ERRORS
	dLIR_dk = np.pi * 4.*np.pi*sqtfR**2 * intgt.quad(lambda nu:2.*h*nu**3 * (nu/nu0)**gam * np.log(nu/nu0)  /(c*c * (-1. + np.exp( h*nu/(kb*Td)) )  ), 0.0, nu0 )[0] / Lav

	


	dLIR_dnu0 = np.pi * 4.*np.pi*sqtfR**2 * intgt.quad(lambda nu:-2.*h*gam*nu*(nu/nu0)**(2. + gam) * nu0 / (c*c * (-1. + np.exp( h*nu/(kb*Td)) )  ), 0.0, nu0 )[0] / Lav

	




	dLIR_dT = np.pi * 4.*np.pi*sqtfR**2 * intg.quad(lambda nu:2.*np.exp(h*nu/(kb*Td)) * h*h * nu**4 * min(1., (nu/nu0)**gam)/( c*c * (-1. + np.exp(h*nu/(kb*Td)) )**2 * kb*Td**2 ), 0.0, numicron*5. )[0] / Lav
	#dLIR_dT = 1.

	#dLIR_dLtot = -4.*np.pi*sqtfR**2 * intgt.quad(lambda nu:min(1., (nu/nu0)**(4.+gam)) * np.pi * Bv(nu, Td), 0.0, numicron*5. )[0] / Lav**2

	#dLIR_dfR2 =  8.*np.pi*sqtfR * intgt.quad(lambda nu:min(1., (nu/nu0)**(4.+gam)) * np.pi * Bv(nu, Td), 0.0, numicron*5. )[0] / Lav

	dLIR_dLtot = -4.*np.pi*sqtfR**2 * intgt.quad(lambda nu:min(1., (nu/nu0)**(gam)) * np.pi * Bv(nu, Td), 0.0, numicron*5. )[0] / Lav**2

	dLIR_dfR2 =  16.*np.pi*sqtfR * intgt.quad(lambda nu:min(1., (nu/nu0)**(gam)) * np.pi * Bv(nu, Td), 0.0, numicron*5. )[0] / Lav

	if (fitFcov):
		delLTot = delLfac * Lav/Lfac  #remember we already updated Lav = Lfac*Lav
	else:
		delLTot = Lav * (1. - 6./10.)
		


	DCF = ma.sqrt( (dLIR_dk*delk)**2 + (dLIR_dnu0*delnu0)**2 + (dLIR_dT*delT)**2  + (dLIR_dLtot*delLTot)**2  + (dLIR_dfR2*delfR2)**2 )
	#DLIR = ma.sqrt(											    (dLIR_dT*delT)**2  + (dLIR_dLtot*delLTot)**2  + (dLIR_dfR2*delfR2)**2 )

	##Error on R(td)
	dRdT = -2.*RR/Td
	dRdL = 0.5*RR/Lav
	dRR = ma.sqrt( (dRdT*delT)**2 + (dRdL * delLTot)**2)
	RTd_print = RR/pc2cm
	dRTd_print = dRR/pc2cm

	## integral over IR BB should be equal to Lbol tot!
	# the 16 is from 4pi aeff^2 *sig_d * 4pi R^2*f and Sigd->1/(pia^2) for tau->1
	# Alternatively take the expression for F_nu, multiply by 4piD^2 and integrate over sphere and freq, if Qnu=1, get LIR/L=f
	Nrm   = 16.*ma.pi*RR**2 * intgt.quad( lambda nu:min(1., (nu/nu0)**(gam)) * np.pi * Bv(nu, Td), 0.0, numicron*5. )[0]
	CF    = 16.*ma.pi*sqtfR**2 * intgt.quad( lambda nu:min(1., (nu/nu0)**(gam)) * np.pi * Bv(nu, Td), 0.0, numicron*5. )[0] / Nrm
	CFmin = 16.*ma.pi*sqtfR**2 * intgt.quad( lambda nu:min(1., (nu/nu0)**(gam)) * np.pi * Bv(nu, Td), 0.0, numicron*5. )[0]  / (Lav*14./10.)
	CFmax = 16.*ma.pi*sqtfR**2 * intgt.quad( lambda nu:min(1., (nu/nu0)**(gam)) * np.pi * Bv(nu, Td), 0.0, numicron*5. )[0]  /(Lav*6./10.)
	CFErrP = CFmax - CF 
	CFErrM = CF - CFmin





	#
	#Rmatch = ma.sqrt(  Lav / (8. * ma.pi  * qIR * h/c/c * (kb/h)**(4.+gam) * spc.gamma(4+gam) * spc.zetac(4.+gam) * Td**(4.+gam)) )
	#Rmatch = sqtfR/ma.sqrt(CF)

	Rmatch = sqtfR/ma.sqrt(Fcov)
	#Rmatch = RR
	Rmatch = Rmatch/pc2cm

	delR = delR/pc2cm


	#delR = ma.sqrt((-0.5*sqtfR/CF**(3./2.) * DCF)**2 + (delfR2/ma.sqrt(CF))**2) / pc2cm

	tdoP = Rmatch*pc2cm/(1474*3600*24 * c)
	delTdoP = tdoP - (Rmatch - delR)*pc2cm/(1474*3600*24 * c)
	nu0_print = nu0/10**14
	delnu0_print = delnu0/10**14
	delfR2 = delfR2/pc2cm
	sqtfR = sqtfR/pc2cm

	nrm = Nrm/Lav

	print "T_opt = %g +- %g K" %(Td, delT)
	print "nu0 = %g +- %g Hz/10^14" %(nu0_print, delnu0_print)
	print "k = %g +- %g" %(gam, delk)

	
	print "sqrt(f)R_d = %g +- %g cm" %(sqtfR, delfR2)

	

	print "R_opt = %g +- %g pc" %(Rmatch, delR)
	print "t_d/P = %g +- %g" %(tdoP, delTdoP)

	print "R_Td = %g +- %g pc" %(RTd_print, dRTd_print)

	print "Int(f=1)/Lav = %g" %nrm


	#TofRopt = 1800.*(  2.*Rmatch/ma.sqrt(CF)/ma.sqrt(Lav/10.**(46)) )**(-1./2.8)
		


	#CF = 0.5*(CFmax + CFmin)
	#CFerr = 0.5*(CFmax - CFmin)

	#print "Fw1_mean = %g +- %g" %(Fw1, ErrW1)

	#print "Fw2_mean = %g +- %g" %(Fw2, ErrW2)

	#print "f = %g +- %g" %(CF, CFerr)
	#print "f = %g + %g - %g" %(CF, CFErrP, CFErrM)
	

	print "Lfac = %g - %g + %g" % (Lfac, dLM, dLP)

	print "f = %g +- %g" %(CF, DCF)

	print "Fcover fit = %g - %g + %g" %(Fcov, dfM, dfP)


	print "chi2 = %g" %chi2_pdf_BB



	### WRITE RESULTS TO FILE
	filename = "../emcee_Results/BBfit_resuts_%iwalkers"%clen+Shell_File+"_FULL_RESULTS.txt" 
	target = open(filename, 'w')
	target.truncate()
	target.write(Shell_File)
	target.write("\n")	
	target.write("T_opt = %g +- %g K" %(Td, delT))
	target.write("\n")	
	target.write("BC = %g - %g + %g" % (Lfac, dLM, dLP))
	target.write("\n")	
	target.write("Fcover fit = %g - %g + %g" %(Fcov, dfM, dfP))
	target.write("\n")	
	target.write("R_opt = %g +- %g pc" %(Rmatch, delR) )
	target.write("\n")	
	target.write( "t_d/P = %g +- %g" %(tdoP, delTdoP) )
	target.write("\n")	
	target.write("chi2 = %g" %chi2_pdf_BB)
	target.close()



	########PLOTTING

	nu = np.linspace(0.01*numicron, 1*numicron, 100)/10**14
	
	pref = np.ones(len(nu))
	for i in range(len(nu)):
		pref[i] = min(1., (nu[i]*10**14/nu0)**(gam))

	plt.figure()
	plt.scatter(W1_mid/10**14, Fw1* 10.**(26), color='orange', s=40, marker='o')
	plt.scatter(W2_mid/10**14, Fw2* 10.**(26), color='red', s =40, marker='o')
	plt.errorbar(W1_mid/10**14, Fw1* 10.**(26), yerr=ErrW1*10.**(26), linestyle="none", color='orange', alpha=1., elinewidth=1.5)
	plt.errorbar(W2_mid/10**14, Fw2* 10.**(26), yerr=ErrW2*10.**(26), linestyle="none", color='red', alpha=1., elinewidth=1.5)

	if (AddW4):
		plt.scatter(W3_mid/10**14, Fw3* 10.**(26), color='purple', s=40, marker='o')
		plt.scatter(W4_mid/10**14, Fw4* 10.**(26), color='brown', s =40, marker='o', alpha=1.0)
		plt.errorbar(W3_mid/10**14, Fw3* 10.**(26), yerr=ErrW3*10.**(26), linestyle="none", color='purple', alpha=1., elinewidth=1.5)
		plt.errorbar(W4_mid/10**14, Fw4* 10.**(26), yerr=ErrW4*10.**(26), linestyle="none", color='brown', alpha=1.0, elinewidth=1.5)
	elif(AddW3):
		plt.scatter(W3_mid/10**14, Fw3* 10.**(26), color='purple', s=40, marker='o')
		plt.scatter(W4_mid/10**14, Fw4* 10.**(26), color='brown', s =40, marker='o', alpha=0.5)
		plt.errorbar(W3_mid/10**14, Fw3* 10.**(26), yerr=ErrW3*10.**(26), linestyle="none", color='purple', alpha=1., elinewidth=1.5)
		plt.errorbar(W4_mid/10**14, Fw4* 10.**(26), yerr=ErrW4*10.**(26), linestyle="none", color='brown', alpha=0.5, elinewidth=1.5)
	else:
		plt.scatter(W3_mid/10**14, Fw3* 10.**(26), color='purple', s=40, marker='o', alpha=0.5)
		plt.scatter(W4_mid/10**14, Fw4* 10.**(26), color='brown', s =40, marker='o', alpha=0.5)
		plt.errorbar(W3_mid/10**14, Fw3* 10.**(26), yerr=ErrW3*10.**(26), linestyle="none", color='purple', alpha=0.5, elinewidth=1.5)
		plt.errorbar(W4_mid/10**14, Fw4* 10.**(26), yerr=ErrW4*10.**(26), linestyle="none", color='brown', alpha=0.5, elinewidth=1.5)


	#plt.axvline(W3_mid/10**14,  color='brown')
	#plt.axvline(W4_mid/10**14,  color='purple')


	#plt.axvline(W1mn/10**14,  color='orange')
	#plt.axvline(W1mx/10**14,  color='orange')
	#plt.axvline(W2mn/10**14,  color='red')
	#plt.axvline(W2mx/10**14,  color='red')



	plt.plot(nu, pref * Bv(nu*10**14, Td) * 4.*ma.pi* (sqtfR*pc2cm/Dst)**2 * 10.**(26), color = 'gray', linewidth = 2)

	plt.xlabel(r'$\nu$ [$10^{14}$ Hz]')
	plt.ylabel('Flux [mJy]')
	
	plt.xlim(0.0,1.5)
	plt.ylim(0.0,110.0)

	plt.tight_layout()

	#plt.show()
	if (fitQv):
		if (AddW3):
			plt.savefig("../emcee_data/BBfit_fitQv_BestFit_ModBlackBody_clen%g_%gwalkersAddW3.png" %(clen, nwalkers))
		else:	
			plt.savefig("../emcee_data/BBfit_fitQv_BestFit_ModBlackBody_clen%g_%gwalkers.png" %(clen, nwalkers))
	else:
		if (AddW4):
			plt.savefig("../emcee_data/BBfit_BestFit_ModBlackBody_clen%g_%gwalkers_AddW4.png" %(clen, nwalkers))
		if (AddW3):
			plt.savefig("../emcee_data/BBfit_BestFit_ModBlackBody_clen%g_%gwalkers_AddW3.png" %(clen, nwalkers))
		else:
			plt.savefig("../emcee_data/BBfit_BestFit_ModBlackBody_clen%g_%gwalkers_W1W2.png" %(clen, nwalkers))


## BB FIT ABOVE

### ANALYTIC PLOTTING BELOW





































if (ISOvDop):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)
	print "Look up tables Created"


	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg1_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]


	Iso_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) 
	Dop_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )

	#LogISO_Anl = np.log10( (1. + Iso_anal) / (1. - Iso_anal) )
	#LogDOP_Anl = np.log10( (1. + Dop_anal) / (1. - Dop_anal) ) 

	ISO_AIR_over_AUV = np.zeros(len(frc_t))
	DOP_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		ISO_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
		DOP_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)

	#Dop_mean = np.mean(DOP_AIR_over_AUV)
	plt.figure()
	AnlIso = plt.plot(frc_a, Iso_anal, color='black', linewidth=3)
	AnlDop = plt.plot(frc_a, Dop_anal, color='red', linewidth=3)
	#AnlIso = plt.plot(frc_a, LogISO_Anl, color='black')
	#AnlDop = plt.plot(frc_a, LogDOP_Anl, color='red')
	ISO = plt.scatter(frc_t, ISO_AIR_over_AUV, marker='x', color='black', s=40)
	DOP = plt.scatter(frc_t, DOP_AIR_over_AUV, marker='*', color='red', s=40)

	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	#plt.axhline(y=Dop_mean, color='red', linestyle='--')

	plt.legend( [ AnlIso[0], AnlDop[0], ISO, DOP ], ('Iso. Analytic', 'Dop. Analytic', 'Iso. Numerical',  'Dop. Numerical'), loc='upper right', fontsize=18)

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)
	plt.ylim(-0.5,1.0)
	#plt.show()

	plt.tight_layout()

	Savename = "plots/Iso_and_Dop/Analytics/DopvsISO_AIRoAUV_alpha%g_TsubCut%g_J%g_numin%g_numx%g_reclim2_TRHS3.png" %(alph, Tsub, JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)




if (ISOvDop_MAX):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin,Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)


	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg1_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]


	Iso_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) 
	Dop_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )

	#LogISO_Anl = np.log10( (1. + Iso_anal) / (1. - Iso_anal) )
	#LogDOP_Anl = np.log10( (1. + Dop_anal) / (1. - Dop_anal) ) 

	ISO_AIR_over_AUV = np.zeros(len(frc_t))
	DOP_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		ISO_AIR_over_AUV[i] = ISO_MAX(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
		DOP_AIR_over_AUV[i] = DOP_MAX(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)

	#Dop_mean = np.mean(DOP_AIR_over_AUV)
	plt.figure()
	plt.title("Obscured IR Minima")
	AnlIso = plt.plot(frc_a, Iso_anal, color='black')
	AnlDop = plt.plot(frc_a, Dop_anal, color='red')
	#AnlIso = plt.plot(frc_a, LogISO_Anl, color='black')
	#AnlDop = plt.plot(frc_a, LogDOP_Anl, color='red')
	ISO = plt.scatter(frc_t, ISO_AIR_over_AUV, marker='x', color='black')
	DOP = plt.scatter(frc_t, DOP_AIR_over_AUV, marker='*', color='red')

	plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	#plt.axhline(y=Dop_mean, color='red', linestyle='--')

	plt.legend( [ AnlIso[0], AnlDop[0], ISO, DOP ], ('Iso Analytic', 'Dop Analytic', 'Iso Max',  'Dop Max'), loc='upper right', fontsize=18)

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)

	plt.tight_layout()

	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/MAXObscDopvsISO_AIRoAUV_J%g_numin%g_numx%g_reclim2_TRHS.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)













if (PG1302_ISO):

	## ranges of plots
	ymx=1.0
	ymn=-1.0
	yc = abs(ymn)/(ymx-ymn)

	A_mx = 0.66
	A_mn = 0.24
	phs_mx = 0.36
	phs_mn = 0.1
	

	if (W1W2vals):
		fmx = 0.05+0.02
		fmn = 0.05-0.02

		tdoPmx = 1.5+0.3
		tdoPmn = 1.5-0.3 

		thT1 = np.arccos(1.)
		
		cosThTmx = fmx
		cosThTmn = fmn
		thT2 = np.arccos(cosThTmn)
		thT3 = np.arccos(cosThTmx)
		JJt = 1.08 #MIN J
		J2  = 1.31
		#J3  = 1.45
		#J4  = 1.4
		J5  = 1.425 #MAX J

		##Find Js where EO Rings bracket allowed region

		##Find Js where Face On Torus at known theta_T bracket allow region

	else:
		fmx = 0.08+0.04
		fmn = 0.08-0.03

		tdoPmx = 3.1+0.8
		tdoPmn = 3.1-0.8 

		#thT1 = 0.0
		thT1 = np.arccos(1)
		cosThTmx = fmx
		cosThTmn = fmn
		thT2 = np.arccos(cosThTmn)
		thT3 = np.arccos(cosThTmx)

		#JJt = 1.435  #MIN J onlypositive AIR/A
		JJt = 1.275
		J2  = 1.435
		#J3  = 1.45
		#J4  = 1.4
		J5  = 1.513 #MAX J


	AA1_anal =  1./(2.*ma.pi*frc_a*ma.cos(thT1)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT1)) 
	AA2_anal =  1./(2.*ma.pi*frc_a*ma.cos(thT2)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT2)) 
	#AA3_anal =  1./(2.*ma.pi*frc_a*ma.cos(thT3)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT3)) 
	AAB_anal =  1./(2.*ma.pi*frc_a*ma.cos(thT3)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT3)) 

	from scipy import special as spc
	Ring_anal = spc.j0(2.*ma.pi*frc_a)


	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)

	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg2_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, J2, aeff, nu0, nne]
	#arg3_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, J3, aeff, nu0, nne]
	#arg4_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, J4, aeff, nu0, nne]
	arg5_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, J5, aeff, nu0, nne]



	Ring_ISO1 = np.zeros(len(frc_t))
	Ring_ISO2 = np.zeros(len(frc_t))
	Ring_ISO3 = np.zeros(len(frc_t))
	Ring_ISO4 = np.zeros(len(frc_t))
	Ring_ISO5 = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		Ring_ISO1[i] = ISO_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
		Ring_ISO2[i] = ISO_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_ISO, RHS_table, T_table)
		#Ring_ISO3[i] = ISO_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_ISO, RHS_table, T_table)
		#Ring_ISO4[i] = ISO_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg4_ISO, RHS_table, T_table)
		Ring_ISO5[i] = ISO_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg5_ISO, RHS_table, T_table)











	plt.figure(figsize=(8,8))

	#plt.title(r'Isotropic, $ \bar{A} / A_{\rm{V}} =%g$, $J = \pi/2$' %Afac)
	#plt.title('Isotropic')
	Anl1   = plt.plot(frc_a, AA1_anal, color='blue', linewidth=3)
	Anl2   = plt.plot(frc_a, AA2_anal, color='blue', linestyle='--', linewidth=1)
	#Anl3   = plt.plot(frc_a, AA3_anal, color='red')
	Ring   = plt.plot(frc_a, Ring_anal, color='green', linewidth=3)
	AnlBst = plt.plot(frc_a, AAB_anal, color='blue', linestyle='--', linewidth=1)

	RingJ1  = plt.plot(frc_t, Ring_ISO1, color='green', linestyle='--',linewidth=1)
	RingJ2   = plt.plot(frc_t, Ring_ISO2, color='green', linestyle='--',linewidth=1)
	#RingJ3   = plt.plot(frc_t, Ring_ISO3, color='green', linewidth=1)
	#RingJ4   = plt.plot(frc_t, Ring_ISO4, color='green', linewidth=1)
	RingJ5   = plt.plot(frc_t, Ring_ISO5, color='green', linestyle='--',linewidth=1)

	#Fillbetween the two limiting Ring solns 
	#plt.fill_between(frc_t, Ring_ISO1, Ring_ISO5, color='green', alpha=0.3)
	plt.fill_between(frc_t, Ring_ISO2, Ring_ISO5, color='green', alpha=0.3)
	#plt.fill_between(frc_t, Ring_ISO1, Ring_ISO2, color='green', alpha=0.3)

	#Fillbetween the two limiting face-on solns 
	plt.fill_between(frc_a, AAB_anal, AA2_anal, color='blue', alpha=0.3)




	##plot td/P measured from L.C. using full W1 and W2 span of phase
	tdoP_W1_mn = phs_mn
	tdoP_W1_mx = phs_mx
	#A_IR>0 
	plt.axvspan(   tdoP_W1_mn, tdoP_W1_mx,    ymin=yc, ymax=1.0, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+1., tdoP_W1_mx+1., ymin=yc, ymax=1.0, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+2., tdoP_W1_mx+2., ymin=yc, ymax=1.0, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+3., tdoP_W1_mx+3., ymin=yc, ymax=1.0, color='orange', alpha=0.4, lw=0, hatch="+")
	#A_IR<0 - 1/2 cycle out
	plt.axvspan(   tdoP_W1_mn + 0.5,    tdoP_W1_mx + 0.5, ymin=0.0, ymax=yc, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+1. + 0.5, tdoP_W1_mx+1. + 0.5, ymin=0.0, ymax=yc, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+2. + 0.5, tdoP_W1_mx+2. + 0.5, ymin=0.0, ymax=yc, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+3. + 0.5, tdoP_W1_mx+3. + 0.5, ymin=0.0, ymax=yc, color='orange', alpha=0.4, lw=0, hatch="+")



	### plot range of td/P from Spectral fit
	plt.axvspan(tdoPmn, tdoPmx, color='orange', alpha=0.5, lw=0)

	### plot amplitude ranges from light curve analysis
	plt.axhspan(A_mn, A_mx, color='red', alpha=0.4, lw=0)
	plt.axhspan(-A_mn, -A_mx, color='red', alpha=0.4, lw=0)

	##plot measured sublimation region
	plt.axvspan(0.0,frc_PG, color='grey', alpha=0.8, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')

	

####################
### OLDER
####################

	###plot td/P Measured from LIR
	#plt.axvline(x=3.821, color='grey', linestyle='--', linewidth=3 )
	#plt.axvline(x=3.3, color='orange', linestyle=':', linewidth=3 )
	#plt.axvspan(3.3-0.7, 3.3+0.7, color='orange', alpha=0.5, lw=0)
	#plt.axvline(x=4.0, color='orange', linestyle=':', linewidth=3 )



	# plot location of Wein peak
	#plt.axvline(x=tdoP_W1_Wein, color='yellow', linestyle=':', linewidth=3 )
	#plt.axvline(x=tdoP_W2_Wein, color='red', linestyle=':', linewidth=3)


	# tdoP_W2_mn = 0.2
	# tdoP_W2_mx = 0.36
	# #A_IR>0 
	# plt.axvspan(tdoP_W2_mn,    tdoP_W2_mx,    ymin=0.333, ymax=1.0, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+1., tdoP_W2_mx+1., ymin=0.333, ymax=1.0, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+2., tdoP_W2_mx+2., ymin=0.333, ymax=1.0, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+3., tdoP_W2_mx+3., ymin=0.333, ymax=1.0, color='red', alpha=0.2, lw=0)

	# #A_IR<0 - 1/2 cycle out
	# plt.axvspan(   tdoP_W2_mn + 0.5,    tdoP_W2_mx + 0.5, ymin=0.0, ymax=0.333, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+1. + 0.5, tdoP_W2_mx+1. + 0.5, ymin=0.0, ymax=0.333, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+2. + 0.5, tdoP_W2_mx+2. + 0.5, ymin=0.0, ymax=0.333, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+3. + 0.5, tdoP_W2_mx+3. + 0.5, ymin=0.0, ymax=0.333, color='red', alpha=0.2, lw=0)



	# ##plot measured AIR/A
	# AW1 = 0.70/Afac
	# AW1_mn = AW1-0.05
	# AW1_mx = AW1+0.05
	# plt.axhspan(AW1_mn, AW1_mx, color='yellow', alpha=0.4, lw=0)
	# AW2 = 0.62/Afac
	# AW2_mn = AW2 - 0.05
	# AW2_mx = AW2 + 0.05
	# plt.axhspan(AW2_mn, AW2_mx, color='red', alpha=0.4, lw=0)

	# # and negatives
	# plt.axhspan(-AW1_mn, -AW1_mx, color='yellow', alpha=0.4, lw=0)
	# plt.axhspan(-AW2_mn, -AW2_mx, color='red', alpha=0.4, lw=0)




	# ##plot measured AIR/A # withbest fit limits on bol correction divided
	# if (W1W2vals):
	# 	Lfac_min = 0.9-0.3
	# 	Lfac_max = 0.9+0.4
	# else:
	# 	Lfac_min = 0.6 - 0.01
	# 	Lfac_max = 0.6 + 0.7



	# AW2 = 0.64/Afac
	# AW2_mn = (0.64-0.20)/2.63/(Lfac_max)
	# AW2_mx = (0.62 +0.12)/(Lfac_min)
	# plt.axhspan(AW2_mn, AW2_mx, color='red', alpha=0.2, lw=0)

	# and negatives
	# plt.axhspan(-AW1_mn, -AW1_mx, color='yellow', alpha=0.2, lw=0)
	# plt.axhspan(-AW2_mn, -AW2_mx, color='red', alpha=0.2, lw=0)

####################
### OLDER
####################
	


	plt.legend( [ Anl1[0], Ring[0], RingJ1[0], AnlBst[0]], (r'$\cos{\theta_T}=1$, $J=\pi/2$', r'$\cos{\theta_T} = 0$' ,r'$\cos{\theta_T} = 0$, $%g \leq J \leq %g$' %(JJt, J5),  r"$%g \leq \cos{\theta_T} \leq %g$, $J=\pi/2$" %(cosThTmn,cosThTmx) ), loc='lower left', fontsize=18)



	plt.ylabel(r"$A_{\rm{IR}} / A$", fontsize=22)
	plt.xlabel(r"$t_d / P$", fontsize=22)



	#plt.xlim(0.0,frc_mx)
	if (W1W2vals):
		plt.xlim(0.0,frc_mx)
	else:
		plt.xlim(0.0,frc_mx)
	plt.ylim(ymn,ymx)

	plt.tight_layout()

	if (W1W2vals):	
		Savename = "plots/Iso_and_Dop/Analytics/PG1302_ISO_W1W2_J%g_numin%g_numx%g.png" %(JJt, Nnumn, Nnumx)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
	else:
		Savename = "plots/Iso_and_Dop/Analytics/PG1302_ISO_W1W2W3_J%g_numin%g_numx%g.png" %(JJt, Nnumn, Nnumx)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)






if (PG1302_Dop):

	Dopalph = -1.0 #Use correct alpha for dust (D^4)
	xmx=1.0
	xmn=-0.7
	yc = abs(xmn)/(xmx-xmn)

	A_mx = 0.38
	A_mn = 0.20
	phs_mx = 0.36
	phs_mn = 0.1
	

	if (W1W2vals):
		fmx = 0.05+0.02
		fmn = 0.05-0.02

		tdoPmx = 1.5+0.3
		tdoPmn = 1.5-0.3 

		thT1 = np.arccos(1)

		#thT3 = np.arccos(0.05+0.02)
		cosThTmn = fmn
		cosThTmx = fmx

		thT2 = np.arccos(cosThTmn)
		thT3 = np.arccos(cosThTmx)
		JJt = 1.43

		#J2  = 1.52
		
		#J3  = 1.52
		#J4  = 1.49
		JRing5 = 1.518
		
		J5  = 1.518

	else:
		fmx = 0.08+0.04
		fmn = 0.08-0.03

		tdoPmx = 3.1+0.8
		tdoPmn = 3.1-0.8 

		thT1 = np.arccos(1)

		#thT3 = np.arccos(0.05+0.02)
		cosThTmn = fmn
		cosThTmx = fmx

		thT2 = np.arccos(cosThTmn)
		thT3 = np.arccos(cosThTmx)
		JJt = 1.49

		#J2  = 1.49
		
		#J3  = 1.52
		#J4  = 1.49
		JRing5 = 1.548
		
		J5  =  1.548



	Inc = 0.0*ma.pi/2.1
	#thT = np.arccos(0.1)


	#JJt = ma.pi/2. - np.arccos(0.1)  #0.1 radians for numerical ring only
	

	
	if (numRing):
		## plot minimum J ring at minimum theta_T
		#nne = 1.8 #no absorption efficiency
		#nu0 = numicron*0.5

		arg1_DOP = [Lav, betst, Inc, Ombn, Dopalph, n0, Rde, pp, thT2, JJt, aeff, nu0, nne]
		#arg2_DOP = [Lav, betst, Inc, Ombn, Dopalph, n0, Rde, pp, thT3, J2, aeff, nu0, nne]
		#arg3_DOP = [Lav, betst, Inc, Ombn, Dopalph, n0, Rde, pp, thT3, J3, aeff, nu0, nne]
		#arg4_DOP = [Lav, betst, Inc, Ombn, Dopalph, n0, Rde, pp, thT3, J4, aeff, nu0, nne]
		arg5_DOP = [Lav, betst, Inc, Ombn, Dopalph, n0, Rde, pp, thT2, J5, aeff, nu0, nne]
		argRing5_DOP = [Lav, betst, Inc, Ombn, Dopalph, n0, Rde, pp, thT2, JRing5, aeff, nu0, nne]

		#argR_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thT3, JJt, aeff, nu0, nne]

		##TABULATE T's and RHSs
		print "Creating look up tables"
		RHS_table = np.zeros(NT)
		T_table = np.linspace(Tmin, Tsub, NT)
		for i in range(NT):
			RHS_table[i] = T_RHS(T_table[i], nu0, nne)



		DOP_AIR_over_AUV = np.zeros(len(frc_t))
		DOP_Ring_AIR_over_AUV = np.zeros(len(frc_t))
		DOP2_Ring_AIR_over_AUV = np.zeros(len(frc_t))
		DOP3_Ring_AIR_over_AUV = np.zeros(len(frc_t))
		DOP4_Ring_AIR_over_AUV = np.zeros(len(frc_t))
		DOP5_Ring_AIR_over_AUV = np.zeros(len(frc_t))
		for i in range(len(frc_t)):
			DOPmnJ_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
			#Assumes theta_T-pi/2
			DOP_Ring_AIR_over_AUV[i] = DOP_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
			#DOP2_Ring_AIR_over_AUV[i] = DOP_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_DOP, RHS_table, T_table)
			#DOP3_Ring_AIR_over_AUV[i] = DOP_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_DOP, RHS_table, T_table)
			#DOP4_Ring_AIR_over_AUV[i] = DOP_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg4_DOP, RHS_table, T_table)
			DOP5_Ring_AIR_over_AUV[i] = DOP_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, argRing5_DOP, RHS_table, T_table)
			DOPmxJ_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg5_DOP, RHS_table, T_table)



	# sin(I-pi/2) = - cos(I)  so add neg in difference
	Dop1_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thT1)) - np.sin(2.*ma.pi*frc_a*np.cos(thT1))/np.cos(thT1) )
	Dop2_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thT2)) - np.sin(2.*ma.pi*frc_a*np.cos(thT2))/np.cos(thT2) )
	Dop3_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thT3)) - np.sin(2.*ma.pi*frc_a*np.cos(thT3))/np.cos(thT3) )
	from scipy import special as spc
	DopRing_anal = spc.j1(2.*ma.pi*frc_a)


	plt.figure(figsize=(8,8))

	
	# sphere ref
	Anl1 = plt.plot(frc_a, Dop1_anal, color='blue', linewidth=3)
	# limits on costhT for face on case
	Anl2 = plt.plot(frc_a, Dop2_anal, color='blue',linestyle='--')
	Anl3 = plt.plot(frc_a, Dop3_anal, color='blue', linestyle='--')
	# ring ref
	Ring = plt.plot(frc_a, DopRing_anal, color='green', linewidth=3)
	if (numRing):
		RingJ = plt.scatter(frc_t, DOP_AIR_over_AUV, marker='*', color='purple')
		#Approx_RingJ = plt.plot(frc_t, DOP_AIR_over_AUV,  color='purple', linestyle='--', linewidth=3)
		RingJ = plt.plot(frc_t, DOP_Ring_AIR_over_AUV,  color='green', linestyle='--', linewidth=1)
		#RingJ2 = plt.plot(frc_t, DOP2_Ring_AIR_over_AUV,  color='green', linestyle='-', linewidth=1)
		#RingJ3 = plt.plot(frc_t, DOP3_Ring_AIR_over_AUV,  color='green', linestyle='-', linewidth=1)
		#RingJ4 = plt.plot(frc_t, DOP4_Ring_AIR_over_AUV,  color='green', linestyle='-', linewidth=1)
		RingJ5 = plt.plot(frc_t, DOP5_Ring_AIR_over_AUV,  color='green', linestyle='--', linewidth=1)

		#plt.title(r'Doppler, $\bar{A} / A_{\rm{V}} =%g$, $I = 0.0$' %Afac)
#	else:
		#plt.title(r'Doppler, $\bar{A} / A_{\rm{V}} =%g$, $J = \pi/2$, $I = 0.0$' %Afac)
		#plt.title('Doppler')

	

	#Fillbetween the two limiting Ring solns 
	plt.fill_between(frc_t, DOP5_Ring_AIR_over_AUV, DOP_Ring_AIR_over_AUV, color='green', alpha=0.3)


	#Fillbetween the two limiting face-on solns 
	plt.fill_between(frc_a, Dop2_anal, Dop3_anal, color='blue', alpha=0.3)



	##plot  td/P measured from Light curves, ambiguous mod one period
	## span W2 and W2 ranges
	# Doppler is a quarter cycle out of phase with UV from get-go
	tdoP_W1_mn = phs_mn + 0.25
	tdoP_W1_mx = phs_mx + 0.25
	#A_IR>0
	plt.axvspan(   tdoP_W1_mn,    tdoP_W1_mx, ymin=yc, ymax=1.0, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+1., tdoP_W1_mx+1., ymin=yc, ymax=1.0, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+2., tdoP_W1_mx+2., ymin=yc, ymax=1.0, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+3., tdoP_W1_mx+3., ymin=yc, ymax=1.0, color='orange', alpha=0.4, lw=0, hatch="+")

	#A_IR<0 - 1/2 cycle out
	plt.axvspan(   tdoP_W1_mn + 0.5,    tdoP_W1_mx + 0.5, ymin=0.0, ymax=yc, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+1. + 0.5, tdoP_W1_mx+1. + 0.5, ymin=0.0, ymax=yc, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+2. + 0.5, tdoP_W1_mx+2. + 0.5, ymin=0.0, ymax=yc, color='orange', alpha=0.4, lw=0, hatch="+")
	plt.axvspan(tdoP_W1_mn+3. + 0.5, tdoP_W1_mx+3. + 0.5, ymin=0.0, ymax=yc, color='orange', alpha=0.4, lw=0, hatch="+")



	### PLOT tdop range fro Spec Fit
	plt.axvspan(tdoPmn, tdoPmx, color='orange', alpha=0.5, lw=0)

	### plot ampl ranges
	plt.axhspan(A_mn, A_mx, color='red', alpha=0.4, lw=0)
	plt.axhspan(-A_mn, -A_mx, color='red', alpha=0.4, lw=0)

	##plot measured sublimation region
	plt.axvspan(0.0,frc_PG, color='grey', alpha=0.8, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')



##################
### OLDER
##################
	## plot td/P Measured from LIR
	#plt.axvline(x=3.821, color='grey', linestyle='--', linewidth=3 )
	#plt.axvline(x=3.3, color='orange', linestyle=':', linewidth=3 )
	#plt.axvspan(3.3-0.7, 3.3+0.7, color='orange', alpha=0.5, lw=0)
	#plt.axvline(x=4.0, color='orange', linestyle=':', linewidth=3 )



	# plot location of Wein peak
	#plt.axvline(x=tdoP_W1_Wein, color='yellow', linestyle=':', linewidth=3 )
	#plt.axvline(x=tdoP_W2_Wein, color='red', linestyle=':', linewidth=3 )

#for W2 (removed do a range for all W1 and W2)
	# tdoP_W2_mn = 0.2 + 0.25
	# tdoP_W2_mx = 0.36 + 0.25
	# #A_IR>0
	# plt.axvspan(   tdoP_W2_mn,    tdoP_W2_mx, ymin=0.333, ymax=1.0, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+1., tdoP_W2_mx+1., ymin=0.333, ymax=1.0, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+2., tdoP_W2_mx+2., ymin=0.333, ymax=1.0, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+3., tdoP_W2_mx+3., ymin=0.333, ymax=1.0, color='red', alpha=0.2, lw=0)

	# #A_IR<0  - 1/2 cycle out
	# plt.axvspan(   tdoP_W2_mn + 0.5,    tdoP_W2_mx + 0.5, ymin=0.0, ymax=0.333, color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+1. + 0.5, tdoP_W2_mx+1. + 0.5, ymin=0.0, ymax=0.333,  color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+2. + 0.5, tdoP_W2_mx+2. + 0.5, ymin=0.0, ymax=0.333,  color='red', alpha=0.2, lw=0)
	# plt.axvspan(tdoP_W2_mn+3. + 0.5, tdoP_W2_mx+3. + 0.5, ymin=0.0, ymax=0.333,  color='red', alpha=0.2, lw=0)


	##plot measured AIR/A
	# AW1 = 0.70/Afac
	# AW1_mn = AW1-0.05
	# AW1_mx = AW1+0.05
	# plt.axhspan(AW1_mn, AW1_mx, color='yellow', alpha=0.4, lw=0)
	# AW2 = 0.62/Afac
	# AW2_mn = AW2 - 0.05
	# AW2_mx = AW2 + 0.05
	# plt.axhspan(AW2_mn, AW2_mx, color='red', alpha=0.4, lw=0)

	##plot measured AIR/A
	# if (W1W2vals):
	# 	Lfac_min = 0.9-0.3
	# 	Lfac_max = 0.9+0.4
	# else:
	# 	Lfac_min = 0.6 - 0.01
	# 	Lfac_max = 0.6 + 0.7

	# AW1_mn = (0.68-0.21)/2.63/(Lfac_max)
	# AW1_mx = (0.68 + 0.12)/(Lfac_min)
	# plt.axhspan(AW1_mn, AW1_mx, color='yellow', alpha=0.2, lw=0)
	# AW2 = 0.64/Afac
	# AW2_mn = (0.64-0.20)/2.63/(Lfac_max)
	# AW2_mx = (0.62 +0.12)/(Lfac_min)
	
	

	# and negatives
	#plt.axhspan(-AW1_mn, -AW1_mx, color='yellow', alpha=0.2, lw=0)
	#plt.axhspan(-AW2_mn, -AW2_mx, color='red', alpha=0.2, lw=0)


	# if (numRing):
	# 	plt.legend( [ Anl1[0],  Anl2[0],  Anl3[0], Ring[0], Dop], (r'$\theta_T = 0$, $J=\pi/2$', r'$\theta_T = \pi/4$, $J=\pi/2$',    r'$\theta_T = \pi/3$, $J=\pi/2$', r'$\theta_T = \pi/2$, $J=0$', r'$\cos{\theta_T} = 0.1$, $J=0.1$'), loc='upper right', fontsize=14)
	# else:
		#plt.legend( [ Anl1[0],  Anl2[0],  Anl3[0], Ring[0]], (r'$\theta_T = 0$, $J=\pi/2$', r'$\theta_T = \pi/4$, $J=\pi/2$',    r'$\theta_T = \pi/3$, $J=\pi/2$', r'$\theta_T = \pi/2$, $J=0$'), loc='upper right', fontsize=14)
		#plt.legend( [ Anl1[0],   Ring[0]], (r'$\theta_T = 0$, $J=\pi/2$', r'$\theta_T = \pi/2$, $J=0$'), loc='upper right', fontsize=14)
	# if (W1W2vals):
	# 	plt.legend( [ Anl1[0], Ring[0], Anl3[0]], (r'$\cos{\theta_T}=1$, $J=\pi/2$', r'$\cos{\theta_T} = 0$, $J=0$', r"$\cos{\theta_T}=%g$, $J=\pi/2$" %cosThT), loc='upper right', fontsize=18)
	# 	if (numRing):
	# 		plt.legend( [ Anl1[0], Ring[0], Anl3[0], RingJ[0]], (r'$\cos{\theta_T}=1$, $J=\pi/2$', r'$\cos{\theta_T} = 0$, $J=0$', r"$\cos{\theta_T}=%g$, $J=\pi/2$" %cosThT, r"$\cos{\theta_{T}}=0$, $J\leq%g$" %JJt, r"$\cos{\theta_{T}}=%g$, $J=%g$" %(cosThT, J5)), loc='upper right', fontsize=18)
	# else:
##################
### OLDER
##################


	if (numRing):
		#plt.legend( [ Anl1[0], Ring[0], Anl3[0], RingJ[0], Approx_RingJ[0]], (r'$\cos{\theta_T}=1$, $J=\pi/2$', r'$\cos{\theta_T} = 0$, $0\leq J \leq \pi/2$', r"$\cos{\theta_T}=%g$, $J=\pi/2$" %cosThT, r"$\cos{\theta_{T}}=0$, $J = %g$" %JJt, r"$\cos{\theta_{T}}=%g$, $J=%g$" %(cosThT, J5) ), loc='upper right', fontsize=18)
		plt.legend( [ Anl1[0], Ring[0], RingJ[0], Anl2[0]], (r'$\cos{\theta_T}=1$, $J=\pi/2$', r'$\cos{\theta_T} = 0$' ,r'$\cos{\theta_T} = 0$, $%g \leq J \leq %g$' %(JJt, J5),  r"$%g \leq \cos{\theta_T} \leq %g$, $J=\pi/2$" %(cosThTmn,cosThTmx) ), loc='lower left', fontsize=18)
	else:
		plt.legend( [ Anl1[0], Ring[0], Anl3[0]], (r'$\cos{\theta_T}=1$, $J=\pi/2$', r'$\cos{\theta_T} = 0$, $J=0$', r"$\cos{\theta_T}=%g$, $J=\pi/2$" %cosThTmx), loc='upper right', fontsize=18)




	plt.ylabel(r"$A_{\rm{IR}} / A$", fontsize=22)
	plt.xlabel(r"$t_d / P$", fontsize=22)


	plt.xlim(0.0,frc_mx)
	plt.ylim(xmn,xmx)

	plt.tight_layout()



	if (W1W2vals):
		Savename = "plots/Iso_and_Dop/Analytics/PG1302_DOP_W1W2_J%g_numin%g_numx%g.png" %(JJt, Nnumn, Nnumx)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
	else:
		Savename = "plots/Iso_and_Dop/Analytics/PG1302_DOP_W1W2W3_J%g_numin%g_numx%g.png" %(JJt, Nnumn, Nnumx)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)



































if (PG1302_ISO_Jth_Conts):
	print "ISO cosT J contours"

	A_mx = 0.66
	A_mn = 0.24
	phs_mx = 0.36 # in units of 2pi td/p /2pi = td/P
	phs_mn = 0.1
	

	if (W1W2vals):
		print "W1W2"
		fmx = 0.05+0.02
		fmn = 0.05-0.02

		tdoPmx = 1.5+0.3
		tdoPmn = 1.5-0.3 

		thTmn = np.arccos(fmn)
		thTmx = np.arccos(fmx)

	else:
		print "W1W2W3"
		fmx = 0.08+0.04
		fmn = 0.08-0.03

		tdoPmx = 3.1+0.8
		tdoPmn = 3.1-0.8 


		thTmn = np.arccos(fmn)
		thTmx = np.arccos(fmx)





	

	# ##TABULATE T's and RHSs
	# print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	# for i in range(NT):
	# 	RHS_table[i] = T_RHS(T_table[i], nu0, nne)




	Nfrc = 12
	Ngrd = 10
	
	
	CosthTs = np.linspace(fmn, fmx, Ngrd)
	thTs = np.arccos(CosthTs)
	JJs = np.linspace(0.0, ma.pi/2., Ngrd)



	AIRoAMax = np.zeros([Ngrd,Ngrd])


	if (W1W2vals):
		Nfrc = Nfrc*2
		AIRoAMax_set1 = np.zeros(Nfrc/2)
		AIRoAMax_set2 = np.zeros(Nfrc/2)
		frc_rng1 = np.linspace(1.2, 1.36, Nfrc)
		frc_rand1 = 1.2 + np.random.rand()*(1.36-1.2)

		frc_rng2 = np.linspace(1.6, 1.8, Nfrc)
		frc_rand2 = 1.6 + np.random.rand()*(1.8-1.6)

		for i in range(Ngrd):
			for j in range(Ngrd):
				k=0
				while (k < (Nfrc/2) ):
					arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thTs[i], JJs[j], aeff, nu0, nne]
			 		#greater than 0
			 		AIRoAMax_set1[k] = ISO_AIR_o_AUV(frc_rng1[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
					# less than zero see if statement below
					AIRoAMax_set2[k] = ISO_AIR_o_AUV(frc_rng2[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table)

					# AIRoAMax_set1[k] = np.abs(ISO_AIR_o_AUV(frc_rand1, numn, numx, Dst, arg1_ISO, RHS_table, T_table))
					# AIRoAMax_set2[k] = np.abs(ISO_AIR_o_AUV(frc_rand2, numn, numx, Dst, arg1_ISO, RHS_table, T_table))

						
					if ( ( (AIRoAMax_set1[k] > 0.0) and (AIRoAMax_set1[k] >= A_mn) and (AIRoAMax_set1[k] <= A_mx) ) or  ( (AIRoAMax_set2[k] < 0.0) and (-AIRoAMax_set2[k] >= A_mn) and (-AIRoAMax_set2[k] <= A_mx) ) ):
						AIRoAMax[j][i] = 1.0
						k = (Nfrc/2) ##EXIT LOOP IF WE HAVE something that works
					else:
						k = k+1
					

					print i, j, k
				# max1 = np.max(AIRoAMax_set1)
				# max2 = np.max(AIRoAMax_set2)
				# min1 = np.min(AIRoAMax_set1)
				# min2 = np.min(AIRoAMax_set2)

				# AIRoAMax[j][i] = np.max( [max1,max2] )
				# AIRoAMin[j][i] = np.min( [min1,min2] )
				# if (AIRoAMax[j][i] < A_mn && AIRoAMin[j][i] > Amx):				
				# 	AIRoAMax[j][i] = 0.0

	else:
		Nfrc = Nfrc*4
		AIRoAMax_set1 = np.zeros(Nfrc/4)
		AIRoAMax_set2 = np.zeros(Nfrc/4)
		AIRoAMax_set3 = np.zeros(Nfrc/4)
		AIRoAMax_set4 = np.zeros(Nfrc/4)

		frc_rng1 = np.linspace(2.3, 2.36 , Nfrc)
		frc_rand1 = 2.3 + np.random.rand()*(2.36-2.3)

		frc_rng2 = np.linspace(2.6, 2.86 , Nfrc)
		frc_rand1 = 2.6 + np.random.rand()*(2.86-2.6)

		frc_rng3 = np.linspace(3.1, 3.36 , Nfrc)
		frc_rand1 = 3.1 + np.random.rand()*(3.36-3.1)

		frc_rng4 = np.linspace(3.6, 3.86 , Nfrc)
		frc_rand1 = 3.6 + np.random.rand()*(3.86-3.6)

		for i in range(Ngrd):
			for j in range(Ngrd):
				k=0
				while (k < (Nfrc/4) ):
					arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thTs[i], JJs[j], aeff, nu0, nne]
			 		AIRoAMax_set1[k] = ISO_AIR_o_AUV(frc_rng1[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
					AIRoAMax_set2[k] = ISO_AIR_o_AUV(frc_rng2[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
					AIRoAMax_set3[k] = ISO_AIR_o_AUV(frc_rng3[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
			 		AIRoAMax_set4[k] = ISO_AIR_o_AUV(frc_rng4[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table)

			 	# 	AIRoAMax_set1[k] = np.abs(ISO_AIR_o_AUV(frc_rand1, numn, numx, Dst, arg1_ISO, RHS_table, T_table))
					# AIRoAMax_set2[k] = np.abs(ISO_AIR_o_AUV(frc_rand2, numn, numx, Dst, arg1_ISO, RHS_table, T_table))
					# AIRoAMax_set3[k] = np.abs(ISO_AIR_o_AUV(frc_rand3, numn, numx, Dst, arg1_ISO, RHS_table, T_table))
			 	# 	AIRoAMax_set4[k] = np.abs(ISO_AIR_o_AUV(frc_rand4, numn, numx, Dst, arg1_ISO, RHS_table, T_table))

					if ( ( (AIRoAMax_set1[k] > 0.0) and (AIRoAMax_set1[k] >= A_mn) and (AIRoAMax_set1[k] <= A_mx) ) or  ( (AIRoAMax_set2[k] < 0.0) and (-AIRoAMax_set2[k] >= A_mn) and (-AIRoAMax_set2[k] <= A_mx) )      or       ( (AIRoAMax_set3[k] > 0.0) and (AIRoAMax_set3[k] >= A_mn) and (AIRoAMax_set3[k] <= A_mx) )     or     ( (AIRoAMax_set4[k] < 0.0) and (-AIRoAMax_set4[k] >= A_mn) and (-AIRoAMax_set4[k] <= A_mx) )):
						AIRoAMax[j][i] = 1.0
						k = (Nfrc/4) ##EXIT LOOP IF WE HAVE something that works
					else:
						k = k+1

					print i, j, k





				# 	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thTs[i], JJs[j], aeff, nu0, nne]
			 # 		AIRoAMax_set1[k] = np.abs(ISO_AIR_o_AUV(frc_rng1[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table))
				# 	AIRoAMax_set2[k] = np.abs(ISO_AIR_o_AUV(frc_rng2[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table))
				# 	AIRoAMax_set3[k] = np.abs(ISO_AIR_o_AUV(frc_rng3[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table))
				# 	AIRoAMax_set4[k] = np.abs(ISO_AIR_o_AUV(frc_rng4[k], numn, numx, Dst, arg1_ISO, RHS_table, T_table))

				# 	print i, j, k
				# max1 = np.max(AIRoAMax_set1)
				# max2 = np.max(AIRoAMax_set2)
				# max3 = np.max(AIRoAMax_set3)
				# max4 = np.max(AIRoAMax_set4)

				# AIRoAMax[j][i] = np.max( [max1,max2,max3,max4] )








	## if this is greater than one, we are a go
	# AMaxratMn = np.log10(AIRoAMax/A_mn - 1.)  ##makes < Amn values nans
	# AMaxratMx = np.log10(A_mx/AIRoAMax - 1.)  ##makes > Amx values nans

	# AMaxrat = (10.**(AMaxratMx)+1. * 10.**(AMaxratMn)+1.)/(10.**(AMaxratMx)+1. ) 

	AMaxrat = AIRoAMax


	plt.figure()
	#plt.title(r"$\left[ A_{\rm{IR}}/A \right]_{\rm{max}}$")
	cnt = plt.contourf(CosthTs, JJs, AMaxrat, 100, cmap = "Greys_r")
	lms = plt.contour(CosthTs, JJs, AMaxrat, cmap = "Greys_r", levels=[0.5])#levels = [A_mn, A_mx])
	#plt.clabel(lms, fmt = '%2.2f', colors = 'k', fontsize=14)	
	#cbar = plt.colorbar(cnt)

	plt.xlabel(r"$\cos{\theta_T}$", fontsize=22)
	plt.ylabel(r"$J$", fontsize=22)


	#plt.xlim(0.0,frc_mx)
	#plt.ylim(xmn,xmx)

	plt.tight_layout()

	#plt.show()
	RdSave = Rde/pc2cm
	if (W1W2vals):
		Savename = "plots/Iso_and_Dop/Analytics/PG1302_ISO_W1W2_TAnyl_JcosT_Contours_tdoPSplit_Rd%gpc_Nfrc%g_Ngrd%g_numin%g_numx%g_Reclim%g_Aerr%g_Rerr%g.png" %(RdSave , Nfrc, Ngrd, Nnumn, Nnumx, reclim, myabs, myrel)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
	else:
		Savename = "plots/Iso_and_Dop/Analytics/PG1302_ISO_W1W2W3__TAnyl_JcosT_Contours_tdoPSplit_Rd%gpc_Nfrc%g_Ngrd%g_numin%g_numx%g_Reclim%g_Aerr%g_Rerr%g.png" %(RdSave , Nfrc, Ngrd, Nnumn, Nnumx, reclim, myabs, myrel)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)















if (PG1302_Dop_Jth_Conts):
	print "DOP cosT J contours"
	Dopalph = -1.0 #Use correct alpha for dust (D^4)
	xmx=1.0
	xmn=-0.7
	yc = abs(xmn)/(xmx-xmn)

	A_mx = 0.38
	A_mn = 0.20
	phs_mx = 0.36
	phs_mn = 0.1
	

	if (W1W2vals):
		print "W1W2"
		fmx = 0.05+0.02
		fmn = 0.05-0.02

		tdoPmx = 1.5+0.3
		tdoPmn = 1.5-0.3 

		thTmn = np.arccos(fmn)
		thTmx = np.arccos(fmx)

	else:
		print "W1W2W3"
		fmx = 0.08+0.04
		fmn = 0.08-0.03

		tdoPmx = 3.1+0.8
		tdoPmn = 3.1-0.8 


		thTmn = np.arccos(fmn)
		thTmx = np.arccos(fmx)



	Inc = 0.0*ma.pi/2.1


	

	# ##TABULATE T's and RHSs
	# print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	# for i in range(NT):
	# 	RHS_table[i] = T_RHS(T_table[i], nu0, nne)




	Nfrc = 12
	Ngrd = 10
	

	CosthTs = np.linspace(fmn, fmx, Ngrd)
	thTs = np.arccos(CosthTs)
	JJs = np.linspace(0.0, ma.pi/2., Ngrd)



	AIRoAMax = np.zeros([Ngrd,Ngrd])


	if (W1W2vals):
		AIRoAMax_set1 = np.zeros(Nfrc)
		frc_rng1 = np.linspace(1.35, 1.61, Nfrc)
		frc_rand1 = 1.35 + np.random.rand()*(1.61-1.35)
		for i in range(Ngrd):
			for j in range(Ngrd):
				k=0
				while (k < (Nfrc) ):
					arg1_DOP = [Lav, betst, Inc, Ombn, Dopalph, n0, Rde, pp, thTs[i], JJs[j], aeff, nu0, nne]
		 		 	AIRoAMax_set1[k] = DOP_AIR_o_AUV(frc_rng1[k], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
			 		#AIRoAMax_set1[k] = np.abs(DOP_AIR_o_AUV(frc_rand1, numn, numx, Dst, arg1_DOP, RHS_table, T_table))

					if ( ((AIRoAMax_set1[k] > 0.0)  and  (AIRoAMax_set1[k] >= A_mn) and (AIRoAMax_set1[k] <= A_mx) ) ):
							AIRoAMax[j][i] = 1.0
							k = (Nfrc) ##EXIT LOOP IF WE HAVE something that works
					else:
						k = k+1

					print i, j, k
	else:
		Nfrc = Nfrc*4
		AIRoAMax_set1 = np.zeros(Nfrc/4)
		AIRoAMax_set2 = np.zeros(Nfrc/4)
		AIRoAMax_set3 = np.zeros(Nfrc/4)
		AIRoAMax_set4 = np.zeros(Nfrc/4)

		frc_rng1 = np.linspace(2.35, 2.61 , Nfrc)
		frc_rand1 = 2.35 + np.random.rand()*(2.61-2.35)

		frc_rng2 = np.linspace(2.85, 3.11 , Nfrc)
		frc_rand2 = 2.85 + np.random.rand()*(3.11-2.85)

		frc_rng3 = np.linspace(3.35, 3.61 , Nfrc)
		frc_rand3 = 3.35 + np.random.rand()*(3.61-3.35)

		frc_rng4 = np.linspace(3.85, 3.90 , Nfrc)
		frc_rand4 = 3.85 + np.random.rand()*(3.9-3.85)

		for i in range(Ngrd):
			for j in range(Ngrd):
				k=0
				while (k < (Nfrc/4) ):
					arg1_DOP = [Lav, betst, Inc, Ombn, Dopalph, n0, Rde, pp, thTs[i], JJs[j], aeff, nu0, nne]
		 		 	AIRoAMax_set1[k] = DOP_AIR_o_AUV(frc_rng1[k], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		 		 	AIRoAMax_set2[k] = DOP_AIR_o_AUV(frc_rng2[k], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		 		 	AIRoAMax_set3[k] = DOP_AIR_o_AUV(frc_rng3[k], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		 		 	AIRoAMax_set4[k] = DOP_AIR_o_AUV(frc_rng4[k], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
						
					# AIRoAMax_set1[k] = np.abs(DOP_AIR_o_AUV(frc_rand1, numn, numx, Dst, arg1_DOP, RHS_table, T_table))
		 		#  	AIRoAMax_set2[k] = np.abs(DOP_AIR_o_AUV(frc_rand2, numn, numx, Dst, arg1_DOP, RHS_table, T_table))
		 		#  	AIRoAMax_set3[k] = np.abs(DOP_AIR_o_AUV(frc_rand3, numn, numx, Dst, arg1_DOP, RHS_table, T_table))
		 		#  	AIRoAMax_set4[k] = np.abs(DOP_AIR_o_AUV(frc_rand4, numn, numx, Dst, arg1_DOP, RHS_table, T_table))
						
							
					if ( ( (AIRoAMax_set1[k] > 0.0) and (AIRoAMax_set1[k] >= A_mn) and (AIRoAMax_set1[k] <= A_mx) ) or  ( (AIRoAMax_set2[k] < 0.0) and (-AIRoAMax_set2[k] >= A_mn) and (-AIRoAMax_set2[k] <= A_mx) )      or       ( (AIRoAMax_set3[k] > 0.0) and (AIRoAMax_set3[k] >= A_mn) and (AIRoAMax_set3[k] <= A_mx) )     or     ( (AIRoAMax_set4[k] < 0.0) and (-AIRoAMax_set4[k] >= A_mn) and (-AIRoAMax_set4[k] <= A_mx) )):
						AIRoAMax[j][i] = 1.0
						k = (Nfrc/4) ##EXIT LOOP IF WE HAVE something that works
					else:
						k = k+1

					print i, j, k








	# Nfrc = 20
	# Ngrd = 20
	# frc_rng = np.arange(tdoPmn, tdoPmx, (tdoPmx-tdoPmn)/Nfrc)
	# AIRoAMax_set = np.zeros(Nfrc)
	# CosthTs = np.arange(fmn, fmx, (fmx-fmn)/Ngrd)
	# thTs = np.arccos(CosthTs)
	# JJs = np.arange(0.0, ma.pi/2, ma.pi/2./Ngrd)

	# AIRoAMax = np.zeros([Ngrd,Ngrd])

	# for i in range(Ngrd):
	# 	for j in range(Ngrd):
	# 		for k in range(Nfrc):
	# 			arg1_DOP = [Lav, betst, Inc, Ombn, Dopalph, n0, Rde, pp, thTs[i], JJs[j], aeff, nu0, nne]
	# 		 	AIRoAMax_set[k] = np.abs(DOP_AIR_o_AUV(frc_rng[k], numn, numx, Dst, arg1_DOP, RHS_table, T_table))
	# 			print i, j, k
	# 		AIRoAMax[j][i] = np.max(AIRoAMax_set)

	# ## if this is greater than one, we are a go
	# # AMaxratMn = np.log10(AIRoAMax/A_mn - 1.)  ##makes < Amn values nans
	# # AMaxratMx = np.log10(A_mx/AIRoAMax - 1.)  ##makes > Amx values nans

	# # AMaxrat = (10.**(AMaxratMx)+1. * 10.**(AMaxratMn)+1.)/(10.**(AMaxratMx)+1. )  ## get back AIRoAMax/A_mn, but with nans

	AMaxrat = AIRoAMax


	plt.figure()
	#plt.title(r"$\log_{10}{\right[\left. \frac{A_{\rm{IR}}{A} \right|_{\rm{max}} / \left. \frac{A_{\rm{IR}}{A} \right|_{\rm{lim}} - 1 \left]}$")
	cnt = plt.contourf(CosthTs, JJs, AMaxrat, 100, cmap = "Greys_r")
	#lms = plt.contour(CosthTs, JJs, AMaxrat, cmap = "gist_heat", levels = [A_mn, A_mx])
	lms = plt.contour(CosthTs, JJs, AMaxrat, cmap = "Greys_r", levels = [0.5])
	#plt.clabel(lms, fmt = '%2.2f', colors = 'k', fontsize=14)	
	#cbar = plt.colorbar(cnt)

	plt.xlabel(r"$\cos{\theta_T}$", fontsize=22)
	plt.ylabel(r"$J$", fontsize=22)


	#plt.xlim(0.0,frc_mx)
	#plt.ylim(xmn,xmx)

	plt.tight_layout()

	#plt.show()
	RdSave = Rde/pc2cm
	if (W1W2vals):
		Savename = "plots/Iso_and_Dop/Analytics/PG1302_DOP_W1W2__TAnyl_JcosT_Contours_Rd%g_Nfrc%g_Ngrd%g_numin%g_numx%g_Reclim%g_Aerr%g_Rerr%g.png" %(RdSave, Nfrc, Ngrd, Nnumn, Nnumx, reclim, myabs, myrel)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
	else:
		Savename = "plots/Iso_and_Dop/Analytics/PG1302_DOP_W1W2W3__TAnyl_JcosT_Contours_Rd%g_Nfrc%g_Ngrd%g_numin%g_numx%g_Reclim%g_Aerr%g_Rerr%g.png" %(RdSave, Nfrc, Ngrd, Nnumn, Nnumx, reclim, myabs, myrel)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)





































if (Dop_alphs):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)

	alph1 = 4.0
	alphA = -1.0
	#alph3 = 2.0
	#alph4 = -2.0

	Inc0 = 0.0
	Inc1 = ma.pi/3.
	Inc2 = ma.pi/2.1
	
	arg1_DOP = [Lav, betst, Inc, Ombn, alph1, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg2_DOP = [Lav, betst, Inc0, Ombn, alphA, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg3_DOP = [Lav, betst, Inc1, Ombn, alphA, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]	
	arg4_DOP = [Lav, betst, Inc2, Ombn, alphA, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]	## set thetTst=0 for this spherical analytica case
	#argA_DOP = [Lav, betst, Inc, Ombn, -1.0, n0, Rde, pp, 0.0, JJt, aeff, nu0, nne]


	#Iso_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) 
	#Dop_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thetTst)*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )
	Dop_anal =  - 1./(2.*ma.pi*frc_a)**2 * (  2.*ma.pi*frc_a*np.cos(thetTst)*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )
	#Dop_tst  =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thetTst)*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )




	
	DOP1_AIR_over_AUV = np.zeros(len(frc_t))
	DOP2_AIR_over_AUV = np.zeros(len(frc_t))
	DOP3_AIR_over_AUV = np.zeros(len(frc_t))
	DOP4_AIR_over_AUV = np.zeros(len(frc_t))
	
	for i in range(len(frc_t)):
		DOP1_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		DOP2_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_DOP, RHS_table, T_table)
		#DOP3_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_DOP, RHS_table, T_table)
		DOP4_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg4_DOP, RHS_table, T_table)
		

	# DOPAnyl_AIR_over_AUV = np.zeros(len(frc_a))
	# for i in range(len(frc_t)):
	# 	DOPAnyl_AIR_over_AUV[i] = DOP_Anyl_AIRoA(frc_a[i], numn, numx, Dst, argA_DOP, RHS_table, T_table)

	DOPAnyl_AIR_over_AUV0 = DOP_Anyl_AIRoA(frc_a, numn, numx, Dst, arg2_DOP, RHS_table, T_table)
	DOPAnyl_AIR_over_AUV1 = DOP_Anyl_AIRoA(frc_a, numn, numx, Dst, arg3_DOP, RHS_table, T_table)
	DOPAnyl_AIR_over_AUV2 = DOP_Anyl_AIRoA(frc_a, numn, numx, Dst, arg4_DOP, RHS_table, T_table)



	#Dop_mean = np.mean(DOP_AIR_over_AUV)
	plt.figure()
	#plt.title(r"$\beta = %g$ $I = %g$" %(betst,Inc) )
	AnlDop = plt.plot(frc_a, Dop_anal, color='black', linewidth=3)
	ADop0 = plt.plot(frc_a, DOPAnyl_AIR_over_AUV0 , color='red', linewidth=2, linestyle="--")
	####ADop1 = plt.plot(frc_a, DOPAnyl_AIR_over_AUV1 , color='blue', linewidth=2, linestyle="-")
	ADop2 = plt.plot(frc_a, DOPAnyl_AIR_over_AUV2 , color='blue', linewidth=2, linestyle="--")


	plt.plot(frc_t, DOP1_AIR_over_AUV, linestyle='--', color='black')
	plt.plot(frc_t, DOP2_AIR_over_AUV, linestyle='--', color='red')
	####plt.plot(frc_t, DOP3_AIR_over_AUV, linestyle='--', color='blue')
	plt.plot(frc_t, DOP4_AIR_over_AUV, linestyle='--', color='blue')

	DOP1 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='black', s=40)
	DOP2 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='red', s=40)
	####DOP3 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='blue', s=40)
	DOP4 = plt.scatter(frc_t, DOP4_AIR_over_AUV, marker='*', color='blue', s=40)

	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	#plt.axhline(y=Dop_mean, color='red', linestyle='--')

	#	plt.legend( [  AnlDop[0], DOP1, DOP2, DOP3, DOP4 ], ('Dop Analytic', r'$\bar{\alpha} = 6.0$', r'$\bar{\alpha} = 4.0$',  r'$\bar{\alpha} = 2.0$', r'$\bar{\alpha} = -2.0$'), loc='upper right', fontsize=18)
	plt.legend( [  AnlDop[0], ADop0[0], ADop2[0]], (r'$\alpha_{\rm{bol}} = 4.0$, Analytic', r'$\alpha_{\rm{bol}} = -1.0$, Aprx $I=0$', r'$\alpha_{\rm{bol}} = -1.0$, Aprx $I=\pi/2$'), loc='upper right', fontsize=18)
	


	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)
	#plt.show()

	plt.tight_layout()

	Savename = "plots/Iso_and_Dop/Analytics/DopvsAlpha_AIRoAUV_J%g_numin%g_numx%g_reclim200_bet%g_TRHS.png" %(JJt, Nnumn, Nnumx, betst)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)














if (Dop_alphs_Rings):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)

	alph1 = 4.0
	alphA = -1.0

	#Approx a Face on Dust ring
	JJt = 0.0#ma.pi/2
	thetRing = ma.pi/2.1

	Inc0 = 0.0
	Inc1 = ma.pi/3.
	Inc2 = ma.pi/2.1

	betst = 0.01
	

	## Plot alph=4 and alph-1 EO ring solution 
	##do for face on and edge on binary, 
	##do for analytic and numerical
	arg1_DOP = [Lav, betst, Inc0, Ombn, alph1, n0, Rde, pp, thetRing, JJt, aeff, nu0, nne]

	arg2_DOP = [Lav, betst, Inc0, Ombn, alphA, n0, Rde, pp, thetRing, JJt, aeff, nu0, nne]
	arg3_DOP = [Lav, betst, Inc2, Ombn, alphA, n0, Rde, pp, thetRing, JJt, aeff, nu0, nne]



	




	## NUMERICAL APPROXS	
	DOP1_AIR_over_AUV = np.zeros(len(frc_t))
	DOP2_AIR_over_AUV = np.zeros(len(frc_t))
	DOP3_AIR_over_AUV = np.zeros(len(frc_t))

	
	for i in range(len(frc_t)):
		DOP1_AIR_over_AUV[i] = DOP_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		DOP2_AIR_over_AUV[i] = DOP_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_DOP, RHS_table, T_table)
		DOP3_AIR_over_AUV[i] = DOP_Ring_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_DOP, RHS_table, T_table)
				

	### Analytic
	Dop4_ANyl =  spc.j1(2.*ma.pi*frc_a)

	DOPAnyl_AIR_over_AUV1 = DOP_EORing3_Anyl_AIRoA(frc_a, numn, numx, Dst, arg2_DOP, RHS_table, T_table)
	DOPAnyl_AIR_over_AUV2 = DOP_EORing3_Anyl_AIRoA(frc_a, numn, numx, Dst, arg3_DOP, RHS_table, T_table)



	#Dop_mean = np.mean(DOP_AIR_over_AUV)
	plt.figure()
	plt.title(r"Doppler Ring Approximation $\beta = %g$" %betst )
	Anl4 = plt.plot(frc_a, Dop4_ANyl, color='black', linewidth=3)
	Anlm1_I0 = plt.plot(frc_a, DOPAnyl_AIR_over_AUV1, color='red', linewidth=3)
	Anlm1_I2 = plt.plot(frc_a, DOPAnyl_AIR_over_AUV2, color='blue', linewidth=3)



	plt.plot(frc_t, DOP1_AIR_over_AUV, linestyle='--', color='black')
	plt.plot(frc_t, DOP2_AIR_over_AUV, linestyle='--', color='red')
	plt.plot(frc_t, DOP3_AIR_over_AUV, linestyle='--', color='blue')
	
	Num4 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='black')
	Numm1_I0 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='red')
	Numm1_I2 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='blue')
	


	
	plt.axhline(y=0, color='black', linestyle=':')
	

	#	plt.legend( [  AnlDop[0], DOP1, DOP2, DOP3, DOP4 ], ('Dop Analytic', r'$\bar{\alpha} = 6.0$', r'$\bar{\alpha} = 4.0$',  r'$\bar{\alpha} = 2.0$', r'$\bar{\alpha} = -2.0$'), loc='upper right', fontsize=18)
	plt.legend( [  Anl4[0], Num4, Anlm1_I0[0], Numm1_I0, Anlm1_I2[0], Numm1_I2 ], (r'$\bar{\alpha} = 4.0$, Analytic', "Num.", r'$\bar{\alpha} = -1.0$, Aprx.', r"Num. $I=0$",r'$\bar{\alpha} = -1.0$, Aprx.', r"Num. $I=\pi/2.1$",), loc='lower right', fontsize=14)


	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0, 1.5)

	plt.tight_layout()

	Savename = "plots/Iso_and_Dop/Analytics/DopvsAlpha_AIRoAUV_J%g_numin%g_numx%g_reclim200_bet%g_TRHS.png" %(JJt, Nnumn, Nnumx, betst)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)






























if (ISOvDop_varyI):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron

	thetTst = ma.pi/4.
	JJt     = 0.*ma.pi/4.
	##TABULATE T's and RHSs
	print "Creating look up tables"

	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)

	Inc1 = 0.0
	Inc2 = ma.pi/4.
	Inc3 = ma.pi/2.1

	arg1_DOP = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg2_DOP = [Lav, betst, Inc2, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg3_DOP = [Lav, betst, Inc3, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]


	
	AA_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thetTst))

	Dop1_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc1-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thetTst)*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )

	Dop2_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc2-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thetTst)*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )

	Dop3_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc3-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thetTst)*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )



	DOP1_AIR_over_AUV = np.zeros(len(frc_t))
	DOP2_AIR_over_AUV = np.zeros(len(frc_t))
	DOP3_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		DOP1_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		DOP2_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_DOP, RHS_table, T_table)
		DOP3_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_DOP, RHS_table, T_table)

	
	plt.figure()
	plt.title(r'Doppler, $J=\pi/2$, $\theta_T = \pi/4$')
	#Anl_ISO = plt.plot(frc_a, AA_anal, color='black')
	#Anl_DOP1 = plt.plot(frc_a, Dop1_anal, color='black')
	#Anl_DOP2 = plt.plot(frc_a, Dop2_anal, color='red')
	#Anl_DOP3 = plt.plot(frc_a, Dop3_anal, color='blue')
	DOP1 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='black')
	DOP2 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='red')
	DOP3 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='blue')

	plt.plot(frc_t, DOP1_AIR_over_AUV, linestyle='--', color='black')
	plt.plot(frc_t, DOP2_AIR_over_AUV, linestyle='--', color='red')
	plt.plot(frc_t, DOP3_AIR_over_AUV, linestyle='--', color='blue')

	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')

	plt.legend( [ DOP1,   DOP2,  DOP3 ], (r"$I=0$", r"$I=\pi/4$", r"$I=\pi/2.1$"), loc='upper right', fontsize=18)

	#plt.legend( [ Anl_DOP1[0], DOP1,  Anl_DOP2[0], DOP2, Anl_DOP3[0], DOP3 ], (r"Analytic $I=0$", r"Numerical", r"Analytic $I=\pi/4$", r"Numerical", r"Analytic $I=\pi/2.1$", r"Numerical"), loc='upper right', fontsize=14)
	#plt.legend( [ DOP1, DOP2, DOP3 ], (r"$I=0$ Dop", r"$I=\pi/4$ Dop", r"$I=\pi/2.1$ Dop"), loc='upper right')



	plt.ylabel(r"$A_{\rm{IR}} / A$")
	#plt.ylabel(r"$A_{\rm{IR}} / (A_{\rm{IR}} + A)$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)

	plt.tight_layout()
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/Dop_varyI_AIRoAUV_alph%g_J%g_numin%g_numx%g_reclim2_TRHS.png" %(alph,JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)












if (ISOvthT):
	thT1 = 0.0
	thT2 = ma.pi/4.
	thT3 = ma.pi/3.
	thT4 = ma.pi/2.05

	#JJt  = 0.*ma.pi/4.

	nne = 0.0 #no absorption efficiency
	##TABULATE T's and RHSs
	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)


	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff, nu0, nne]
	arg2_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT2, JJt, aeff, nu0, nne]
	arg3_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT3, JJt, aeff, nu0, nne]
	arg4_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT4, JJt, aeff, nu0, nne]



	AA1_anal =  1./(2.*ma.pi*frc_a* ma.cos(thT1)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT1)) 
	AA2_anal =  1./(2.*ma.pi*frc_a* ma.cos(thT2)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT2)) 
	AA3_anal =  1./(2.*ma.pi*frc_a* ma.cos(thT3)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT3)) 
	AA4_anal =  1./(2.*ma.pi*frc_a* ma.cos(thT4)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT4)) 



	from scipy import special as spc
	Ring_anal = spc.j0(2.*ma.pi*frc_a)

	ISO1_AIR_over_AUV = np.zeros(len(frc_t))
	ISO2_AIR_over_AUV = np.zeros(len(frc_t))
	ISO3_AIR_over_AUV = np.zeros(len(frc_t))
	ISO4_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		ISO1_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
		ISO2_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_ISO, RHS_table, T_table)
		ISO3_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_ISO, RHS_table, T_table)
		ISO4_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg4_ISO, RHS_table, T_table)



	plt.figure()

	if (JJt==ma.pi/2):
		plt.title(r'Isotropic, $J = \pi/2$')
		Anl1 = plt.plot(frc_a, AA1_anal, color='black',  linewidth=3)
		Anl2 = plt.plot(frc_a, AA2_anal, color='blue',  linewidth=3)
		Anl3 = plt.plot(frc_a, AA3_anal, color='red',  linewidth=3)
		Anl4 = plt.plot(frc_a, AA4_anal, color='green',  linewidth=3)

		

		ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black', s=40)
		ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='blue', s=40)
		ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='red', s=40)
		ISO4 = plt.scatter(frc_t, ISO4_AIR_over_AUV, marker='x', color='green', s=40)
		plt.plot(frc_t, ISO1_AIR_over_AUV, linestyle='--', color='black')
		plt.plot(frc_t, ISO2_AIR_over_AUV, linestyle='--', color='blue')
		plt.plot(frc_t, ISO3_AIR_over_AUV, linestyle='--', color='red')
		plt.plot(frc_t, ISO4_AIR_over_AUV, linestyle='--', color='green')




		#plt.axvline(x=frc_PG, color='black', linestyle=':')
		#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
		plt.axhline(y=0, color='black', linestyle=':')


		
		#plt.legend( [ Anl1[0],  ISO1, Anl2[0], ISO2, Anl3[0], ISO3, Anl4[0], ISO4], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$', '',    r'$\theta_T = \pi/3$', '', r'$\theta_T = \pi/2.05$', ''), loc='upper right', fontsize=18)
		plt.legend( [ Anl1[0], Anl2[0], Anl3[0], Anl4[0]], (r'$\theta_T = 0$',  r'$\theta_T = \pi/4$',    r'$\theta_T = \pi/3$',  r'$\theta_T = \pi/2.05$'), loc='upper right', fontsize=18)


	elif (JJt==0.0):
		plt.title(r'Isotropic, $J = 0.0$')
		Anl1 = plt.plot(frc_a, AA1_anal, color='black',  linewidth=3)

		Ring = plt.plot(frc_a, Ring_anal, color='green',  linewidth=3)

		ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black', s=40)
		ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='blue', s=40)
		ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='red', s=40)
		plt.plot(frc_t, ISO1_AIR_over_AUV, linestyle='--', color='black')
		plt.plot(frc_t, ISO2_AIR_over_AUV, linestyle='--', color='blue')
		plt.plot(frc_t, ISO3_AIR_over_AUV, linestyle='--', color='red')


		plt.axhline(y=0, color='black', linestyle=':')


		plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3, Ring[0] ], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$',    r'$\theta_T = \pi/3$', 'Ring'), loc='upper right', fontsize=18)
	

	elif (JJt==ma.pi/4.):
		plt.title(r'Isotropic, $J = \pi/4$')
		Anl1 = plt.plot(frc_a, AA1_anal, color='black',  linewidth=3)

		ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black', s=40)
		ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='blue', s=40)
		ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='red', s=40)
		ISO4 = plt.scatter(frc_t, ISO4_AIR_over_AUV, marker='x', color='green', s=40)
		plt.plot(frc_t, ISO1_AIR_over_AUV, linestyle='--', color='black')
		plt.plot(frc_t, ISO2_AIR_over_AUV, linestyle='--', color='blue')
		plt.plot(frc_t, ISO3_AIR_over_AUV, linestyle='--', color='red')
		plt.plot(frc_t, ISO4_AIR_over_AUV, linestyle='--', color='green')


		plt.axhline(y=0, color='black', linestyle=':')


		plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3, ISO4 ], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$',    r'$\theta_T = \pi/3$',    r'$\theta_T = \pi/2.05$'), loc='upper right', fontsize=18)

	else:
		print "Pick another J ISO"




	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")



	plt.xlim(0.0,frc_mx)
	plt.ylim(-0.5,1.0)

	plt.tight_layout()
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/AIRoAUV_ylimmp5_REL_J%g_numin%g_numx%g_reclim1.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)








if (DOPvthT):
	thT1 = 0.0
	thT2 = ma.pi/4.
	thT3 = ma.pi/3.
	thT4 = ma.pi/2.05

	#JJt  = 0.*ma.pi/2.

	Inc = 0.0*ma.pi/2.1

	nne = 0.0 #no absorption efficiency
	##TABULATE T's and RHSs
	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)

	arg1_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thT1, JJt, aeff, nu0, nne]
	arg2_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thT2, JJt, aeff, nu0, nne]
	arg3_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thT3, JJt, aeff, nu0, nne]
	arg4_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thT4, JJt, aeff, nu0, nne]



	Dop1_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thT1)) - np.sin(2.*ma.pi*frc_a*np.cos(thT1))/np.cos(thT1) )

	Dop2_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thT2)) - np.sin(2.*ma.pi*frc_a*np.cos(thT2))/np.cos(thT2) )

	Dop3_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thT3)) - np.sin(2.*ma.pi*frc_a*np.cos(thT3))/np.cos(thT3) )
	
	Dop4_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thT4)) - np.sin(2.*ma.pi*frc_a*np.cos(thT4))/np.cos(thT4) )


	from scipy import special as spc
	DopRing_anal = spc.j1(2.*ma.pi*frc_a)

	DOP1_AIR_over_AUV = np.zeros(len(frc_t))
	DOP2_AIR_over_AUV = np.zeros(len(frc_t))
	DOP3_AIR_over_AUV = np.zeros(len(frc_t))
	DOP4_AIR_over_AUV = np.zeros(len(frc_t))


	for i in range(len(frc_t)):
		DOP1_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		DOP2_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_DOP, RHS_table, T_table)
		DOP3_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_DOP, RHS_table, T_table)
		DOP4_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg4_DOP, RHS_table, T_table)


	plt.figure()

	if (JJt == ma.pi/2.):
		plt.title(r'Doppler, $J = \pi/2$, $I = %g$' %Inc)


		Anl1 = plt.plot(frc_a, Dop1_anal, color='black',  linewidth=3)
		Anl2 = plt.plot(frc_a, Dop2_anal, color='blue',  linewidth=3)
		Anl3 = plt.plot(frc_a, Dop3_anal, color='red',  linewidth=3)
		Anl4 = plt.plot(frc_a, Dop4_anal, color='green',  linewidth=3)
		DOP1 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='black', s=40)
		DOP2 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='blue', s=40)
		DOP3 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='red', s=40)
		DOP4 = plt.scatter(frc_t, DOP4_AIR_over_AUV, marker='*', color='green', s=40)
		plt.plot(frc_t, DOP1_AIR_over_AUV, color='black', linestyle="--")
		plt.plot(frc_t, DOP2_AIR_over_AUV, color='blue', linestyle="--")
		plt.plot(frc_t, DOP3_AIR_over_AUV, color='red', linestyle="--")
		plt.plot(frc_t, DOP4_AIR_over_AUV, color='green', linestyle="--")


		#plt.axvline(x=frc_PG, color='black', linestyle=':')
		#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
		plt.axhline(y=0, color='black', linestyle=':')

		#plt.legend( [ Anl1[0],  DOP1, Anl2[0], DOP2, Anl3[0], DOP3, Anl4[0], DOP4 ], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$', '',   r'$\theta_T = \pi/3$', '',   r'$\theta_T = \pi/2.05$', ''), loc='upper right', fontsize=18)

		plt.legend( [ Anl1[0],   Anl2[0],  Anl3[0],  Anl4[0] ], (r'$\theta_T = 0$',   r'$\theta_T = \pi/4$',    r'$\theta_T = \pi/3$',    r'$\theta_T = \pi/2.05$'), loc='upper right', fontsize=18)
			
	elif (JJt == 0.0):
		plt.title(r'Doppler, $J = 0.0$, $I = %g$' %Inc)


		Anl1 = plt.plot(frc_a, Dop1_anal, color='black',  linewidth=3)
		Ring = plt.plot(frc_a, DopRing_anal, color='green',  linewidth=3)


		DOP1 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='black', s=40)
		DOP2 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='blue', s=40)
		DOP3 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='red', s=40)
		plt.plot(frc_t, DOP1_AIR_over_AUV, color='black', linestyle="--")
		plt.plot(frc_t, DOP2_AIR_over_AUV, color='blue', linestyle="--")
		plt.plot(frc_t, DOP3_AIR_over_AUV, color='red', linestyle="--")


		
		plt.axhline(y=0, color='black', linestyle=':')


		
		plt.legend( [ Anl1[0],  DOP1,  DOP2,  DOP3 , Ring[0]], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$',  r'$\theta_T = \pi/3$', 'Ring'), loc='upper right', fontsize=18)
	elif (JJt == ma.pi/4.):
		plt.title(r'Doppler, $J = \pi/4$, $I = %g$' %Inc)


		Anl1 = plt.plot(frc_a, Dop1_anal, color='black',  linewidth=3)

		DOP1 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='black',s=40)
		DOP2 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='blue', s=40)
		DOP3 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='red', s=40)
		DOP4 = plt.scatter(frc_t, DOP4_AIR_over_AUV, marker='*', color='green', s=40)
		plt.plot(frc_t, DOP1_AIR_over_AUV, color='black', linestyle="--")
		plt.plot(frc_t, DOP2_AIR_over_AUV, color='blue', linestyle="--")
		plt.plot(frc_t, DOP3_AIR_over_AUV, color='red', linestyle="--")
		plt.plot(frc_t, DOP4_AIR_over_AUV, color='green', linestyle="--")

		plt.axhline(y=0, color='black', linestyle=':')


		plt.legend( [ Anl1[0],  DOP1,  DOP2,  DOP3, DOP4 ], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$',  r'$\theta_T = \pi/3$',  r'$\theta_T = \pi/2.05$'), loc='upper right', fontsize=18)



		plt.ylabel(r"$A_{\rm{IR}} / A$")
		plt.xlabel(r"$t_d / P$")



		plt.xlim(0.0,frc_mx)
		plt.ylim(-0.5,1.0)

		plt.tight_layout()
		#plt.show()

		Savename = "plots/Iso_and_Dop/Analytics/AIRoAUVDop_REL_J%g_numin%g_numx%g_reclim2.png" %(JJt, Nnumn, Nnumx)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
















if (Qv_k):
	thT1 = 0.0
	JJt  = ma.pi/2.

	nn1 = 0.0
	print "Creating look up tables"
	RHS_table1 = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table1[i] = T_RHS(T_table[i], nu0, nn1)


	nn2 = 1.0 #no absorption efficiency
	print "Creating look up tables"
	RHS_table2 = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table2[i] = T_RHS(T_table[i], nu0, nn2)

	nn3 = 2.0 #no absorption efficiency
	print "Creating look up tables"
	NT = 15000
	RHS_table3 = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table3[i] = T_RHS(T_table[i], nu0, nn3)

	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff, nu0, nn1]
	arg2_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff, nu0, nn2]
	arg3_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff, nu0, nn3]



	AA1_anal =  1./(2.*ma.pi*frc_a* ma.cos(thT1)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT1)) 
	#AA2_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT2)) 
	#AA3_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT3)) 

	ISO1_AIR_over_AUV = np.zeros(len(frc_t))
	ISO2_AIR_over_AUV = np.zeros(len(frc_t))
	ISO3_AIR_over_AUV = np.zeros(len(frc_t))


	for i in range(len(frc_t)):
		ISO1_AIR_over_AUV[i]  = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table1, T_table)
		ISO2_AIR_over_AUV[i]  = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_ISO, RHS_table2, T_table)
		ISO3_AIR_over_AUV[i]  = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_ISO, RHS_table3, T_table)



	plt.figure()
	plt.title(r'Isotropic')
	#plt.title(r'$J = \pi/2$')
	Anl1 = plt.plot(frc_a, AA1_anal, color='black',  linewidth=3)

	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black', s=40)
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='red', s=40)
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='blue', s=40)
	
	

	#plt.axvline(x=frc_PG, color='black', linestyle=':')
	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	

	plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3], ("Analytic", r"$k=0$", r"$k=1$", r"$k=2$"), loc='upper right', fontsize=18)

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	

	plt.xlim(0.0,frc_mx)

	plt.tight_layout()
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/Qv_vark__AIRoAUV_J%g_numin%g_numx%g_reclim1.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)




if (Qv_nu0):
	thT1 = 0.0
	JJt  = ma.pi/2.




	#kk = 10.0
	kk = 1.6
	fac1 = 1./2.8
	nu01 = fac1 * numicron
	aeff1 = (c/nu01)/(2.*ma.pi)

	fac2 = 2./(2.8 + 4.)
	nu02 = fac2 * numicron
	aeff2 = (c/nu02)/(2.*ma.pi)

	fac3 = 1./4.
	nu03 = fac3 * numicron
	aeff3 = (c/nu03)/(2.*ma.pi)

	Tsub = 3000.
	NT = 10000

	print "Creating nu01 look up tables"
	RHS_table1 = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table1[i] = T_RHS(T_table[i], nu01, kk)


	
	print "Creating nu02 look up tables"
	RHS_table2 = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table2[i] = T_RHS(T_table[i], nu02, kk)

	
	print "Creating nu03 look up tables"
	RHS_table3 = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table3[i] = T_RHS(T_table[i], nu03, kk)


# #### get approx peak BB freq at average
# arg1_T = [Lav, 0.0, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff1, numicron, 10]


# kb = 1.3807*10**(-16)
# h = 6.62607*10**(-27) 
# nupk1 = 2.821439 * kb/h * TDust_Iso(0., Rde, ma.pi/2., 0.0, arg1_T, RHS_table1, T_table)
	







	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff1, nu01, kk]
	arg2_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff2, nu02, kk]
	arg3_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff3, nu03, kk]



	AA1_anal =  1./(2.*ma.pi*frc_a* ma.cos(thT1)) * np.sin(2.*ma.pi*frc_a * ma.cos(thT1)) 
	

	ISO1_AIR_over_AUV = np.zeros(len(frc_t))
	ISO2_AIR_over_AUV = np.zeros(len(frc_t))
	ISO3_AIR_over_AUV = np.zeros(len(frc_t))


	for i in range(len(frc_t)):
		ISO1_AIR_over_AUV[i]  = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table1, T_table)
		ISO2_AIR_over_AUV[i]  = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_ISO, RHS_table2, T_table)
		ISO3_AIR_over_AUV[i]  = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_ISO, RHS_table3, T_table)



	plt.figure()
	plt.title(r'Isotropic')
	#plt.title(r'$J = \pi/2$')
	Anl1 = plt.plot(frc_a, AA1_anal, color='black',  linewidth=3)

	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black', s=40)
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='red', s=40)
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='blue', s=40)



	#plt.axvline(x=frc_PG, color='black', linestyle=':')
	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')


	#plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3], ("Analytic", r"$%g\nu_{\mu m}$" %fac1, r"$%g\nu_{\mu m}$" %fac2, r"$%g\nu_{\mu m}$" %fac3), loc='upper right', fontsize=18)
	plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3], ("Analytic", r"$2.8 \mu m$", r"$3.4 \mu m$", r"$4.0 \mu m$"), loc='upper right', fontsize=18)


	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")



	plt.xlim(0.0,frc_mx)

	plt.tight_layout()
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/Qv_varnu0_SublimationCutT%g_AIRoAUV_k%g_J%g_numin%g_numx%g_reclim1.png" %(Tsub, kk, JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)











if (ThkvThn_ISO):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)


	
	
	JJ1 = 0.0
	JJ2 = ma.pi/4.
	JJ3 = ma.pi/2

	thetT1 = ma.pi/4.
	#thetT2 = ma.pi/4.
	#thetT3 = ma.pi/4.

	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetT1, JJ1, aeff, nu0, nne]
	arg2_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetT1, JJ2, aeff, nu0, nne]
	arg3_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetT1, JJ3, aeff, nu0, nne]



	AA_anal =  1./(2.*ma.pi*frc_a* ma.cos(thetT1)) * np.sin(2.*ma.pi*frc_a * ma.cos(thetT1)) 

	ISO1_AIR_over_AUV = np.zeros(len(frc_t))
	ISO2_AIR_over_AUV = np.zeros(len(frc_t))
	ISO3_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		ISO1_AIR_over_AUV[i] = ISO_OptThick_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
		ISO2_AIR_over_AUV[i] = ISO_OptThick_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_ISO, RHS_table, T_table)
		ISO3_AIR_over_AUV[i] = ISO_OptThick_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_ISO, RHS_table, T_table)

		#DOP_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)

	#Dop_mean = np.mean(DOP_AIR_over_AUV)
	plt.figure()
	plt.title(r'Isotropic, $\theta_T = \pi/4$')
	Anl = plt.plot(frc_a, AA_anal, color='black')
	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black')
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='red')
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='blue')
	plt.plot(frc_t, ISO1_AIR_over_AUV, marker='x', color='black')
	plt.plot(frc_t, ISO2_AIR_over_AUV, marker='x', color='red')
	plt.plot(frc_t, ISO3_AIR_over_AUV, marker='x', color='blue')


	#DOP = plt.scatter(frc_t, DOP_AIR_over_AUV, marker='*', color='red')

	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	#plt.axhline(y=Dop_mean, color='red', linestyle='--')

	#plt.legend( [ Anl[0], ISO, DOP ], ('Iso Analytic', 'Iso numerical',  'Dop numerical'), loc='upper right')

	plt.legend( [ Anl[0], ISO1, ISO2, ISO3 ], ('Iso Opt. Thin', r'Opt.Thick $J=0$', r'$J=\pi/4$', r'$J=\pi/2$'), loc='upper right', fontsize=18)


	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)
	plt.ylim(-0.5,1.0)

	Savename = "plots/Iso_and_Dop/Analytics/Thick_vThin_ISO_AIRoAUV_J%g_numin%g_numx%g_reclim2_TRHS.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)







if (ThkvThn_Dop):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"

	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)


	
	
	JJ1 = 0.0
	JJ2 = ma.pi/4.
	JJ3 = ma.pi/2

	Inc = 0.0

	thetT1 = ma.pi/4.
	#thetT2 = ma.pi/4.
	#thetT3 = ma.pi/4.

	arg1_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetT1, JJ1, aeff, nu0, nne]
	arg2_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetT1, JJ2, aeff, nu0, nne]
	arg3_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetT1, JJ3, aeff, nu0, nne]



	Dop_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst))/np.cos(thetTst) )

	DOP1_AIR_over_AUV = np.zeros(len(frc_t))
	DOP2_AIR_over_AUV = np.zeros(len(frc_t))
	DOP3_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		DOP1_AIR_over_AUV[i] = DOP_OptThick_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		DOP2_AIR_over_AUV[i] = DOP_OptThick_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_DOP, RHS_table, T_table)
		DOP3_AIR_over_AUV[i] = DOP_OptThick_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_DOP, RHS_table, T_table)

		#DOP_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)

	#Dop_mean = np.mean(DOP_AIR_over_AUV)
	plt.figure()
	plt.title(r'Doppler, $\theta_T = \pi/4$')
	Anl = plt.plot(frc_a, Dop_anal, color='black')
	DOP1 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='black')
	DOP2 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='red')
	DOP3 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='blue')
	plt.plot(frc_t, DOP1_AIR_over_AUV, marker='*', color='black')
	plt.plot(frc_t, DOP2_AIR_over_AUV, marker='*', color='red')
	plt.plot(frc_t, DOP3_AIR_over_AUV, marker='*', color='blue')

	#DOP = plt.scatter(frc_t, DOP_AIR_over_AUV, marker='*', color='red')

	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	#plt.axhline(y=Dop_mean, color='red', linestyle='--')

	#plt.legend( [ Anl[0], ISO, DOP ], ('Iso Analytic', 'Iso numerical',  'Dop numerical'), loc='upper right')

	plt.legend( [ Anl[0], DOP1, DOP2, DOP3 ], ('Dop Opt. Thin', r'Opt.Thick $J=0$', r'$J=\pi/4$', r'$J=\pi/2$'), loc='upper right', fontsize=18)


	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)
	plt.ylim(-0.5,1.0)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/Thick_vThin_DOP_AIRoAUV_J%g_numin%g_numx%g_reclim2_TRHS.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)
























































if (fit_Dop):

	from scipy import optimize
	from scipy.optimize import fmin

	def Err2(p, RdOPc, y, dy):
		print p
		AA, fac = p
		chi = (y -  AA*np.sin(fac*RdOPc*2.*ma.pi)/(fac*RdOPc*2.*ma.pi) ) / dy
		chi2 = sum(chi*chi)
		print chi2
		return chi2

	print "Creating look up tables"

	RHS_table = np.zeros(NT)
	T_table = np.linspace(Tmin, Tsub, NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, 0.0)


	arg_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	print "Filling Amplitude array"
	RdOPcs = np.linspace(0.01, 2.0, 20.0)
	ys = np.zeros(len(RdOPcs))
	for i in range(0, len(RdOPcs)):
		ys[i] = DOP_AIR_o_AUV(RdOPcs[i], numn, numx, Dst, arg_DOP, RHS_table, T_table)

	# ys =[ 0.00843292,  0.12321739,  0.18851526,  0.08520175, -0.07485579,
 #       -0.1192214 , -0.06182829,  0.01618426,  0.01383309, -0.04385433,
 #       -0.05974861, -0.03515941,  0.00329961,  0.00314962, -0.0252608 ,
 #       -0.03961657, -0.0232662 ,  0.00042709, -0.00758721, -0.02597535]

	dys = 0.0001*RdOPcs/RdOPcs


	print "Fitting"
	p0 = [ 0.11483526,  7.11715477]
	p_opt = [ 0.11483526,  7.11715477]
	#p_opt  = sc.optimize.fmin(Err2,    p0, args=(RdOPcs, ys, dys), full_output=1, disp=False,ftol=0.000001)[0]


	print p_opt




	print "Plotting"
	plt.figure()

	

	Dop = plt.scatter(RdOPcs, ys, marker='x', color='black')
	fit = plt.plot(RdOPcs, p_opt[0]*np.sin(p_opt[1]*RdOPcs)/(p_opt[1]*RdOPcs), color='black')
	plt.axhline(y=0, color='black', linestyle='--')
	

	plt.legend( [ fit[0],  Dop], ("Fit", "Numerical"), loc='upper right')

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	

	plt.xlim(0.0,2.0)
	plt.show()


