import matplotlib
#matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] = 'sans-serif'
#matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 16})


from FluxFuncs_IRLE import *






###OPTIONS###OPTIONS
ISOvDop = False
ISOvDop_MAX = False

Dop_alphs = False
ISOvDop_varyI = False

ISOvthT = False
PG1302_ISO = True
PG1302_Dop = True

DOPvthT = False

Qv_nu0 = False
Qv_k = False


ThkvThn_ISO = False
ThkvThn_Dop = False

fit_Dop = False





##ISO VARS:
Amp = 0.22105  ##For alpha=0.0, beta=0.0677
t0 = ma.pi/2. #0.0




#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
Omb = 1./(1*yr2sec)
L0 = 6.78*10**46
MPGmx = 10**9.4*Msun
Ryr = c*yr2sec
RdPG = ma.sqrt(0.1)*2.8 *pc2cm
OmPG = 2.*ma.pi/(1474*3600*24) #Omb*2.*ma.pi/4.1
Ombn = OmPG
#alphnu = 1.1

Rorb = c*2.*ma.pi/Omb
Ompc = 2.*ma.pi*c/pc2cm/2.


## TEST VALUES
### DUST stuff
## for Qv
nne = 1.
nu0 = numicron

Rde = RdPG
Rrout = 1.0*Rde
pp = 2.0
thetTst = 0.0#1.*ma.pi/2.
JJt = ma.pi/2.

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




## INTEGRATION LIMTS FOR ALL nu
#Nnumn = 1./2.8 #0.0#0.00001
#Nnumx = 1./4.0 #5.0

Nnumn = 0.0
Nnumx = 5.0

numn = Nnumn*numicron
numx = Nnumx*numicron




frc_PG = Rde /ma.sqrt(0.1)  * ma.sqrt( 10**(8.7)/(10**(9.0)) )/c / (2.*ma.pi/OmPG) ## sublimation radisu for L = Lav (10^8.7 is Mass for which epsilon = 1)
frc_PGb = frc_PG /ma.sqrt(0.1)

frc_mx = 3.0
frc_t = np.linspace(0.01, frc_mx, 30.0)
frc_a = np.linspace(0.01, frc_mx, 100.0)



## First Compute amplitdue
def ISO_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc
	t0 = 0.0
	arg2 = [Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	LIR_mx = F_Sphere_Iso_QuadInt(numn, numx, 0.25*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Sphere_Iso_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff
	

	#LIR_bol = F_Sphere_Iso_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	#return np.abs(0.5*(LIR_mx - LIR_mn)/(Lav*Amp))
	
	return 0.5*(LIR_mx - LIR_mn)/(Lav*Amp) #*   Lav/(0.5*(LIR_mx + LIR_mn))
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
	#nne = 0.0
	#thetTst = 
	arg2 = [Lav, betst, Inc, Ombn, alph, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	
	LUV_mx = Fsrc_Dop(0.5*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LUV_mn = Fsrc_Dop(0.75*2*ma.pi/Ombn, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
	LIR_mx = F_Sphere_Dop_QuadInt(numn, numx, 0.5*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Sphere_Dop_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	

	#LIR_bol = F_Sphere_Dop_QuadInt(0.0, 5.*numicron, 0.5*2*ma.pi/Ombn, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 


	#return np.abs((LIR_mx - LIR_mn)/(LUV_mx - LUV_mn))
	return (LIR_mx - LIR_mn)/(LUV_mx - LUV_mn) #* Lav/(0.5*(LIR_mx + LIR_mn))
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

	return (LIR_mx - LIR_mn)/(LUV_mx - LUV_mn) 


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












if (ISOvDop):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	Tsub = 2000.
	NT = 10000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., Tsub, NT)
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
		ISO_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
		DOP_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)

	#Dop_mean = np.mean(DOP_AIR_over_AUV)
	plt.figure()
	AnlIso = plt.plot(frc_a, Iso_anal, color='black')
	AnlDop = plt.plot(frc_a, Dop_anal, color='red')
	#AnlIso = plt.plot(frc_a, LogISO_Anl, color='black')
	#AnlDop = plt.plot(frc_a, LogDOP_Anl, color='red')
	ISO = plt.scatter(frc_t, ISO_AIR_over_AUV, marker='x', color='black')
	DOP = plt.scatter(frc_t, DOP_AIR_over_AUV, marker='*', color='red')

	plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	#plt.axhline(y=Dop_mean, color='red', linestyle='--')

	plt.legend( [ AnlIso[0], AnlDop[0], ISO, DOP ], ('Iso. Analytic', 'Dop. Analytic', 'Iso. Numerical',  'Dop. Numerical'), loc='upper right', fontsize=14)

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)
	#plt.ylim(-0.5,1.0)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/DopvsISO_AIRoAUV_alpha%g_TsubCut%g_J%g_numin%g_numx%g_reclim50_TRHS2mx.png" %(alph, Tsub, JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)




if (ISOvDop_MAX):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	NT = 15000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 3000., NT)
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

	plt.legend( [ AnlIso[0], AnlDop[0], ISO, DOP ], ('Iso Analytic', 'Dop Analytic', 'Iso Max',  'Dop Max'), loc='upper right', fontsize=14)

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)

	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/MAXObscDopvsISO_AIRoAUV_J%g_numin%g_numx%g_reclim2_TRHS.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)













if (PG1302_ISO):

	Afac = 1.0
	thT1 = 0.0
	thT2 = ma.pi/4.
	thT3 = ma.pi/3.


	AA1_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT1)) 
	AA2_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT2)) 
	AA3_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT3)) 

	

	plt.title(r'Isotropic, $ \bar{A} / A_{\rm{V}} =%g$, $J = \pi/2$' %Afac)
	Anl1 = plt.plot(frc_a, AA1_anal, color='black')
	Anl2 = plt.plot(frc_a, AA2_anal, color='blue')
	Anl3 = plt.plot(frc_a, AA3_anal, color='red')




	##plot measured sublimation region
	#plt.axvline(x=frc_PG, color='black', linestyle=':')

	plt.axvspan(frc_PG, frc_mx, color='grey', alpha=0.4, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')


	##plot measured td/P
	tdoP_W1_mn = 0.1
	tdoP_W1_mx = 0.29
	plt.axvspan(tdoP_W1_mn, tdoP_W1_mx, color='yellow', alpha=0.4, lw=0)
	plt.axvspan(tdoP_W1_mn+1., tdoP_W1_mx+1., color='yellow', alpha=0.4, lw=0)
	plt.axvspan(tdoP_W1_mn+2., tdoP_W1_mx+2., color='yellow', alpha=0.4, lw=0)

	tdoP_W2_mn = 0.2
	tdoP_W2_mx = 0.36
	plt.axvspan(tdoP_W2_mn, tdoP_W2_mx, color='red', alpha=0.4, lw=0)
	plt.axvspan(tdoP_W2_mn+1., tdoP_W2_mx+1., color='red', alpha=0.4, lw=0)
	plt.axvspan(tdoP_W2_mn+2., tdoP_W2_mx+2., color='red', alpha=0.4, lw=0)

	##plot measured AIR/A
	AW1 = 0.70/Afac
	AW1_mn = AW1-0.05
	AW1_mx = AW1+0.05
	plt.axhspan(AW1_mn, AW1_mx, color='yellow', alpha=0.4, lw=0)
	AW2 = 0.62/Afac
	AW2_mn = AW2 - 0.05
	AW2_mx = AW2 + 0.05
	plt.axhspan(AW2_mn, AW2_mx, color='red', alpha=0.4, lw=0)


	
	plt.legend( [ Anl1[0],  Anl2[0],  Anl3[0] ], (r'$\theta_T = 0$', r'$\theta_T = \pi/4$',    r'$\theta_T = \pi/3$'), loc='upper right', fontsize=14)



	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")



	plt.xlim(0.0,frc_mx)
	plt.ylim(-0.4,1.0)


	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/PG1302_ISO_divAfac%g_J%g_numin%g_numx%g_reclim1.png" %(Afac, JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)






if (PG1302_Dop):

	Afac = 1.0
	thT1 = 0.0
	thT2 = ma.pi/4.
	thT3 = ma.pi/3.


	Inc = 0.0*ma.pi/2.1

	
	Dop1_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thT1)*np.cos(2.*ma.pi*frc_a*np.cos(thT1)) - np.sin(2.*ma.pi*frc_a*np.cos(thT1)) )
	Dop2_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thT2)*np.cos(2.*ma.pi*frc_a*np.cos(thT2)) - np.sin(2.*ma.pi*frc_a*np.cos(thT2)) )
	Dop3_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thT3)*np.cos(2.*ma.pi*frc_a*np.cos(thT3)) - np.sin(2.*ma.pi*frc_a*np.cos(thT3)) )


	plt.figure()

	plt.title(r'Doppler, $ \bar{A} / A_{\rm{V}} =%g$, $J = \pi/2$, $I = 0.0$' %Afac)

	Anl1 = plt.plot(frc_a, Dop1_anal, color='black')
	Anl2 = plt.plot(frc_a, Dop2_anal, color='blue')
	Anl3 = plt.plot(frc_a, Dop3_anal, color='red')


	


	
	##plot measured sublimation region
	plt.axvspan(frc_PG, frc_mx, color='grey', alpha=0.4, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')


	##plot measured td/P
	tdoP_W1_mn = 0.1
	tdoP_W1_mx = 0.29
	plt.axvspan(tdoP_W1_mn, tdoP_W1_mx, color='yellow', alpha=0.4, lw=0)
	plt.axvspan(tdoP_W1_mn+1., tdoP_W1_mx+1., color='yellow', alpha=0.4, lw=0)
	plt.axvspan(tdoP_W1_mn+2., tdoP_W1_mx+2., color='yellow', alpha=0.4, lw=0)

	tdoP_W2_mn = 0.2
	tdoP_W2_mx = 0.36
	plt.axvspan(tdoP_W2_mn, tdoP_W2_mx, color='red', alpha=0.4, lw=0)
	plt.axvspan(tdoP_W2_mn+1., tdoP_W2_mx+1., color='red', alpha=0.4, lw=0)
	plt.axvspan(tdoP_W2_mn+2., tdoP_W2_mx+2., color='red', alpha=0.4, lw=0)

	##plot measured AIR/A
	AW1 = 0.70/Afac
	AW1_mn = AW1-0.05
	AW1_mx = AW1+0.05
	plt.axhspan(AW1_mn, AW1_mx, color='yellow', alpha=0.4, lw=0)
	AW2 = 0.62/Afac
	AW2_mn = AW2 - 0.05
	AW2_mx = AW2 + 0.05
	plt.axhspan(AW2_mn, AW2_mx, color='red', alpha=0.4, lw=0)


	
	plt.legend( [ Anl1[0],  Anl2[0],  Anl3[0] ], (r'$\theta_T = 0$', r'$\theta_T = \pi/4$',    r'$\theta_T = \pi/3$'), loc='upper right', fontsize=14)



	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")



	plt.xlim(0.0,frc_mx)
	plt.ylim(-0.4,1.0)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/PG1302_DOP_divAfac%g_J%g_numin%g_numx%g_reclim2.png" %(Afac, JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)






































if (Dop_alphs):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	NT = 15000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 3000., NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)

	alph1 = 6.0
	alph2 = 4.0
	alph3 = 2.0
	alph4 = -2.0

	Inc = 0.0
	
	arg1_DOP = [Lav, betst, Inc, Ombn, alph1, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg2_DOP = [Lav, betst, Inc, Ombn, alph2, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg3_DOP = [Lav, betst, Inc, Ombn, alph3, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg4_DOP = [Lav, betst, Inc, Ombn, alph4, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]


	#Iso_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) 
	Dop_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thetTst)*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )

	
	DOP1_AIR_over_AUV = np.zeros(len(frc_t))
	DOP2_AIR_over_AUV = np.zeros(len(frc_t))
	DOP3_AIR_over_AUV = np.zeros(len(frc_t))
	DOP4_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		DOP1_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		DOP2_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_DOP, RHS_table, T_table)
		DOP3_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_DOP, RHS_table, T_table)
		DOP4_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg4_DOP, RHS_table, T_table)


	#Dop_mean = np.mean(DOP_AIR_over_AUV)
	plt.figure()
	AnlDop = plt.plot(frc_a, Dop_anal, color='red')

	plt.plot(frc_t, DOP1_AIR_over_AUV, linestyle='--', color='black')
	plt.plot(frc_t, DOP2_AIR_over_AUV, linestyle='--', color='red')
	plt.plot(frc_t, DOP3_AIR_over_AUV, linestyle='--', color='blue')
	plt.plot(frc_t, DOP4_AIR_over_AUV, linestyle='--', color='green')

	DOP1 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='black')
	DOP2 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='red')
	DOP3 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='blue')
	DOP4 = plt.scatter(frc_t, DOP4_AIR_over_AUV, marker='*', color='green')

	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	#plt.axhline(y=Dop_mean, color='red', linestyle='--')

	plt.legend( [  AnlDop[0], DOP1, DOP2, DOP3, DOP4 ], ('Dop Analytic', r'$\alpha = 6.0$', r'$\alpha = 4.0$',  r'$\alpha = 2.0$', r'$\alpha = -2.0$'), loc='upper right', fontsize=14)

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/DopvsAlpha_AIRoAUV_J%g_numin%g_numx%g_reclim1_TRHS.png" %(JJt, Nnumn, Nnumx)
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
	NT = 15000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 3000., NT)
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

	plt.legend( [ DOP1,   DOP2,  DOP3 ], (r"$I=0$", r"$I=\pi/4$", r"$I=\pi/2.1$"), loc='upper right', fontsize=14)

	#plt.legend( [ Anl_DOP1[0], DOP1,  Anl_DOP2[0], DOP2, Anl_DOP3[0], DOP3 ], (r"Analytic $I=0$", r"Numerical", r"Analytic $I=\pi/4$", r"Numerical", r"Analytic $I=\pi/2.1$", r"Numerical"), loc='upper right', fontsize=14)
	#plt.legend( [ DOP1, DOP2, DOP3 ], (r"$I=0$ Dop", r"$I=\pi/4$ Dop", r"$I=\pi/2.1$ Dop"), loc='upper right')



	plt.ylabel(r"$A_{\rm{IR}} / A$")
	#plt.ylabel(r"$A_{\rm{IR}} / (A_{\rm{IR}} + A)$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/Dop_varyI_AIRoAUV_alph%g_J%g_numin%g_numx%g_reclim2_TRHS.png" %(alph,JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)












if (ISOvthT):
	thT1 = 0.0
	thT2 = ma.pi/4.
	thT3 = ma.pi/3.
	JJt  = 0.*ma.pi/4.

	nne = 0.0 #no absorption efficiency
	##TABULATE T's and RHSs
	print "Creating look up tables"
	NT = 10000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)


	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff, nu0, nne]
	arg2_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT2, JJt, aeff, nu0, nne]
	arg3_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT3, JJt, aeff, nu0, nne]



	AA1_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT1)) 
	AA2_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT2)) 
	AA3_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT3)) 

	ISO1_AIR_over_AUV = np.zeros(len(frc_t))
	ISO2_AIR_over_AUV = np.zeros(len(frc_t))
	ISO3_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		ISO1_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
		ISO2_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_ISO, RHS_table, T_table)
		ISO3_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_ISO, RHS_table, T_table)



	plt.figure()

	plt.title(r'Isotropic, $J = 0.0$')
	Anl1 = plt.plot(frc_a, AA1_anal, color='black')
	#Anl2 = plt.plot(frc_a, AA2_anal, color='blue')
	#Anl3 = plt.plot(frc_a, AA3_anal, color='red')

	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black')
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='blue')
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='red')
	plt.plot(frc_t, ISO1_AIR_over_AUV, linestyle='--', color='black')
	plt.plot(frc_t, ISO2_AIR_over_AUV, linestyle='--', color='blue')
	plt.plot(frc_t, ISO3_AIR_over_AUV, linestyle='--', color='red')



	#plt.axvline(x=frc_PG, color='black', linestyle=':')
	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')


	#plt.legend( [ Anl1[0],  ISO1, Anl2[0], ISO2, Anl3[0], ISO3 ], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$', '',    r'$\theta_T = \pi/3$', '',), loc='upper right', fontsize=14)
	plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3 ], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$',    r'$\theta_T = \pi/3$'), loc='upper right', fontsize=14)



	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")



	plt.xlim(0.0,frc_mx)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/AIRoAUV_J%g_numin%g_numx%g_reclim1.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)








if (DOPvthT):
	thT1 = 0.0
	thT2 = ma.pi/4.
	thT3 = ma.pi/3.

	JJt  = 0.*ma.pi/4.

	Inc = 0.0*ma.pi/2.1

	nne = 0.0 #no absorption efficiency
	##TABULATE T's and RHSs
	print "Creating look up tables"
	NT = 15000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 3000., NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)

	arg1_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thT1, JJt, aeff, nu0, nne]
	arg2_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thT2, JJt, aeff, nu0, nne]
	arg3_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thT3, JJt, aeff, nu0, nne]


	#AA1_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT1)) 
	#AA2_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT2)) 
	#AA3_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT3)) 


	Dop1_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thT1)*np.cos(2.*ma.pi*frc_a*np.cos(thT1)) - np.sin(2.*ma.pi*frc_a*np.cos(thT1)) )

	Dop2_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thT2)*np.cos(2.*ma.pi*frc_a*np.cos(thT2)) - np.sin(2.*ma.pi*frc_a*np.cos(thT2)) )

	Dop3_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(thT3)*np.cos(2.*ma.pi*frc_a*np.cos(thT3)) - np.sin(2.*ma.pi*frc_a*np.cos(thT3)) )



	DOP1_AIR_over_AUV = np.zeros(len(frc_t))
	DOP2_AIR_over_AUV = np.zeros(len(frc_t))
	DOP3_AIR_over_AUV = np.zeros(len(frc_t))


	for i in range(len(frc_t)):
		DOP1_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)
		DOP2_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg2_DOP, RHS_table, T_table)
		DOP3_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg3_DOP, RHS_table, T_table)



	plt.figure()

	plt.title(r'Doppler, $J = 0.0$, $I = 0.0$')
	#Anl1 = plt.plot(frc_a, AA1_anal, color='black')
	#Anl2 = plt.plot(frc_a, AA2_anal, color='blue')
	#Anl3 = plt.plot(frc_a, AA3_anal, color='red')

	Anl1 = plt.plot(frc_a, Dop1_anal, color='black')
	#Anl2 = plt.plot(frc_a, Dop2_anal, color='blue')
	#Anl3 = plt.plot(frc_a, Dop3_anal, color='red')
	DOP1 = plt.scatter(frc_t, DOP1_AIR_over_AUV, marker='*', color='black')
	DOP2 = plt.scatter(frc_t, DOP2_AIR_over_AUV, marker='*', color='blue')
	DOP3 = plt.scatter(frc_t, DOP3_AIR_over_AUV, marker='*', color='red')
	plt.plot(frc_t, DOP1_AIR_over_AUV, color='black', linestyle="--")
	plt.plot(frc_t, DOP2_AIR_over_AUV, color='blue', linestyle="--")
	plt.plot(frc_t, DOP3_AIR_over_AUV, color='red', linestyle="--")


	#plt.axvline(x=frc_PG, color='black', linestyle=':')
	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')


	#plt.legend( [ Anl1[0],  DOP1, Anl2[0], DOP2, Anl3[0], DOP3 ], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$', '',   r'$\theta_T = \pi/3$', ''), loc='upper right', fontsize=14)
	plt.legend( [ Anl1[0],  DOP1,  DOP2,  DOP3 ], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$',  r'$\theta_T = \pi/3$'), loc='upper right', fontsize=14)



	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")



	plt.xlim(0.0,frc_mx)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/AIRoAUVDop_J%g_numin%g_numx%g_reclim2.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)
















if (Qv_k):
	thT1 = 0.0
	JJt  = ma.pi/2.

	nn1 = 0.0
	print "Creating look up tables"
	NT = 15000
	RHS_table1 = np.zeros(NT)
	T_table = np.linspace(1., 3000., NT)
	for i in range(NT):
		RHS_table1[i] = T_RHS(T_table[i], nu0, nn1)


	nn2 = 1.0 #no absorption efficiency
	print "Creating look up tables"
	NT = 15000
	RHS_table2 = np.zeros(NT)
	T_table = np.linspace(1., 3000., NT)
	for i in range(NT):
		RHS_table2[i] = T_RHS(T_table[i], nu0, nn2)

	nn3 = 2.0 #no absorption efficiency
	print "Creating look up tables"
	NT = 15000
	RHS_table3 = np.zeros(NT)
	T_table = np.linspace(1., 3000., NT)
	for i in range(NT):
		RHS_table3[i] = T_RHS(T_table[i], nu0, nn3)

	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff, nu0, nn1]
	arg2_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff, nu0, nn2]
	arg3_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thT1, JJt, aeff, nu0, nn3]



	AA1_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT1)) 
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
	Anl1 = plt.plot(frc_a, AA1_anal, color='black')

	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black')
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='red')
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='blue')
	
	

	#plt.axvline(x=frc_PG, color='black', linestyle=':')
	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	

	plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3], ("Analytic", r"$k=0$", r"$k=1$", r"$k=2$"), loc='upper right', fontsize=14)

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	

	plt.xlim(0.0,frc_mx)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/Qv_vark__AIRoAUV_J%g_numin%g_numx%g_reclim1.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)




if (Qv_nu0):
	thT1 = 0.0
	JJt  = ma.pi/2.




	kk = 2.0
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
	T_table = np.linspace(1., Tsub, NT)
	for i in range(NT):
		RHS_table1[i] = T_RHS(T_table[i], nu01, kk)


	
	print "Creating nu02 look up tables"
	RHS_table2 = np.zeros(NT)
	T_table = np.linspace(1., Tsub, NT)
	for i in range(NT):
		RHS_table2[i] = T_RHS(T_table[i], nu02, kk)

	
	print "Creating nu03 look up tables"
	RHS_table3 = np.zeros(NT)
	T_table = np.linspace(1., Tsub, NT)
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



	AA1_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a * ma.cos(thT1)) 
	

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
	Anl1 = plt.plot(frc_a, AA1_anal, color='black')

	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black')
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='red')
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='blue')



	#plt.axvline(x=frc_PG, color='black', linestyle=':')
	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')


	plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3], ("Analytic", r"$%g\nu_0$" %fac1, r"$%g\nu_0$" %fac2, r"$%g\nu_0$" %fac3), loc='upper right', fontsize=14)

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")



	plt.xlim(0.0,frc_mx)
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
	NT = 10000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 3000., NT)
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



	AA_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a) 

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

	#DOP = plt.scatter(frc_t, DOP_AIR_over_AUV, marker='*', color='red')

	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	#plt.axhline(y=Dop_mean, color='red', linestyle='--')

	#plt.legend( [ Anl[0], ISO, DOP ], ('Iso Analytic', 'Iso numerical',  'Dop numerical'), loc='upper right')

	plt.legend( [ Anl[0], ISO1, ISO2, ISO3 ], ('Iso Opt. Thin', r'Opt.Thick $J=0$', r'$J=\pi/4$', r'$J=\pi/2$'), loc='upper right', fontsize=14)


	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/Thick_vThin_ISO_AIRoAUV_J%g_numin%g_numx%g_reclim2_TRHS.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)







if (ThkvThn_Dop):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	NT = 10000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 3000., NT)
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



	Dop_anal =  1./(2.*ma.pi*frc_a)**2 *np.sin(Inc-np.pi/2) * (  2.*ma.pi*frc_a*np.cos(2.*ma.pi*frc_a*np.cos(thetTst)) - np.sin(2.*ma.pi*frc_a*np.cos(thetTst)) )

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

	#DOP = plt.scatter(frc_t, DOP_AIR_over_AUV, marker='*', color='red')

	#plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle=':')
	#plt.axhline(y=Dop_mean, color='red', linestyle='--')

	#plt.legend( [ Anl[0], ISO, DOP ], ('Iso Analytic', 'Iso numerical',  'Dop numerical'), loc='upper right')

	plt.legend( [ Anl[0], DOP1, DOP2, DOP3 ], ('Dop Opt. Thin', r'Opt.Thick $J=0$', r'$J=\pi/4$', r'$J=\pi/2$'), loc='upper right', fontsize=14)


	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,frc_mx)
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
	NT = 1000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
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


