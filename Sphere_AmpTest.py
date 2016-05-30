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
ISOvthT = False
Qv_nu0 = True
Qv_k = False


# I_name = "Incs_Snd"

# R_name = "Rdust_Snd"
# Plot_Om  = True
# Om_name = "Ombins_Snd"
# Plot_bet = True
# bet_name = "betas_Snd"
# Plot_J = True
# J_name = "Js_Snd"
# Plot_TT = True
# TT_name = "Theta_Ts_Snd"
# Plot_pp = True
# pp_name = "_ps_Snd"
# Plot_Ro = True
# Ro_name = "_Routs_Snd"
# Plot_n0 = False




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
alph = 0.0
Dst = 1.4*10**9*pc2cm




### WISE BAND + Observational STUFF
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3

nuVbnd = c/(545.*10**(-7))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-20)*(W1mn + W1mx)/2
FW2Rel = 1.7187*10**(-20)*(W2mn + W2mx)/2


# if (not (Qv_k or Qv_nu0)):
# 	##TABULATE T's and RHSs
# 	print "Creating look up tables"
# 	NT = 10000
# 	RHS_table = np.zeros(NT)
# 	T_table = np.linspace(1., 2000., NT)
# 	for i in range(NT):
# 		RHS_table[i] = T_RHS(T_table[i], nu0, nne)




## INTEGRATION LIMTS FOR ALL nu
Nnumn = 0.0#0.00001
Nnumx = 100.0
numn = Nnumn*numicron
numx = Nnumx*numicron

## First Compute amplitdue
def ISO_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
	Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne = arg1
	Ombn = 2. *ma.pi * c/Rd * frc
	t0 = 0.0
	arg2 = [Lav, Amp, Ombn, t0, n0, Rd, pp, thetTst, JJt, aeff, nu0, nne]
	LIR_mx = F_Sphere_Iso_QuadInt(numn, numx, 0.25*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	LIR_mn = F_Sphere_Iso_QuadInt(numn, numx, 0.75*2*ma.pi/Ombn + Rd/c, Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff 
	
	#return np.abs(0.5*(LIR_mx - LIR_mn)/(Lav*Amp))
	return 0.5*(LIR_mx - LIR_mn)/(Lav*Amp)

def DOP_AIR_o_AUV(frc, numn, numx, Dst, arg1, RHS_table, T_table):
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
	return (LIR_mx - LIR_mn)/(LUV_mx - LUV_mn)



frc_PG = Rde/c / (2.*ma.pi/OmPG)
frc_PGb = frc_PG /ma.sqrt(0.1)

frc_t = np.linspace(0.01, 3.0, 30.0)
frc_a = np.linspace(0.01, 3.0, 100.0)


if (ISOvDop):
	nne = 0.0 #no absorption efficiency
	nu0 = numicron
	##TABULATE T's and RHSs
	print "Creating look up tables"
	NT = 10000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)


	arg1_ISO = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
	arg1_DOP = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]


	AA_anal =  1./(2.*ma.pi*frc_a) * np.sin(2.*ma.pi*frc_a) 

	ISO_AIR_over_AUV = np.zeros(len(frc_t))
	DOP_AIR_over_AUV = np.zeros(len(frc_t))
	for i in range(len(frc_t)):
		ISO_AIR_over_AUV[i] = ISO_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_ISO, RHS_table, T_table)
		DOP_AIR_over_AUV[i] = DOP_AIR_o_AUV(frc_t[i], numn, numx, Dst, arg1_DOP, RHS_table, T_table)

	Dop_mean = np.mean(DOP_AIR_over_AUV)
	plt.figure()
	Anl = plt.plot(frc_a, AA_anal, color='black')
	ISO = plt.scatter(frc_t, ISO_AIR_over_AUV, marker='x', color='black')
	DOP = plt.scatter(frc_t, DOP_AIR_over_AUV, marker='*', color='red')

	plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle='--')
	plt.axhline(y=Dop_mean, color='red', linestyle='--')

	plt.legend( [ Anl[0], ISO, DOP ], ('Analytic', 'Iso numerical',  'Dop numerical'), loc='upper right')

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	plt.xlim(0.0,3.0)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/DopvsISO_AIRoAUV_J%g_numin%g_numx%g_reclim2.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)




if (ISOvthT):
	thT1 = 0.0
	thT2 = ma.pi/4.
	thT3 = ma.pi/3.
	JJt  = ma.pi/2.

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

	plt.title(r'$J = 0$')
	Anl1 = plt.plot(frc_a, AA1_anal, color='black')
	Anl2 = plt.plot(frc_a, AA2_anal, color='blue')
	Anl3 = plt.plot(frc_a, AA3_anal, color='red')
	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black')
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='blue')
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='red')
	

	#plt.axvline(x=frc_PG, color='black', linestyle=':')
	plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle='--')
	

	plt.legend( [ Anl1[0],  ISO1, Anl2[0], ISO2, Anl3[0], ISO3 ], (r'$\theta_T = 0$', '',  r'$\theta_T = \pi/4$', '',    r'$\theta_T = \pi/3$', '',), loc='upper right')

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	

	plt.xlim(0.0,3.0)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/AIRoAUV_J%g_numin%g_numx%g_reclim2.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)


if (Qv_k):
	thT1 = 0.0
	JJt  = ma.pi/2.

	nn1 = 0.0
	print "Creating look up tables"
	NT = 10000
	RHS_table1 = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
	for i in range(NT):
		RHS_table1[i] = T_RHS(T_table[i], nu0, nn1)


	nn2 = 1.0 #no absorption efficiency
	print "Creating look up tables"
	NT = 10000
	RHS_table2 = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
	for i in range(NT):
		RHS_table2[i] = T_RHS(T_table[i], nu0, nn2)

	nn3 = 2.0 #no absorption efficiency
	print "Creating look up tables"
	NT = 10000
	RHS_table3 = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
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

	plt.title(r'$J = \pi/2$')
	Anl1 = plt.plot(frc_a, AA1_anal, color='black')

	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black')
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='red')
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='blue')
	
	

	#plt.axvline(x=frc_PG, color='black', linestyle=':')
	plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle='--')
	

	plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3], ("Analytic", r"$k=0$", r"$k=1$", r"$k=2$"), loc='upper right')

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	

	plt.xlim(0.0,3.0)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/Qv_vark__AIRoAUV_J%g_numin%g_numx%g_reclim2.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)




if (Qv_nu0):
	thT1 = 0.0
	thT2 = ma.pi/4.
	thT3 = ma.pi/3.
	JJt  = ma.pi/2.

	kk = 1.0
	nu01 = 0.5 * numicron
	aeff1 = (c/nu01)/(2.*ma.pi)

	nu02 = numicron
	aeff2 = (c/nu02)/(2.*ma.pi)

	nu03 = 2.0 * numicron
	aeff3 = (c/nu01)/(2.*ma.pi)


	print "Creating look up tables"
	NT = 10000
	RHS_table1 = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
	for i in range(NT):
		RHS_table1[i] = T_RHS(T_table[i], nu01, kk)


	
	print "Creating look up tables"
	NT = 10000
	RHS_table2 = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
	for i in range(NT):
		RHS_table2[i] = T_RHS(T_table[i], nu02, kk)

	
	print "Creating look up tables"
	NT = 10000
	RHS_table3 = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
	for i in range(NT):
		RHS_table3[i] = T_RHS(T_table[i], nu03, kk)

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

	plt.title(r'$J = \pi/2$')
	Anl1 = plt.plot(frc_a, AA1_anal, color='black')

	ISO1 = plt.scatter(frc_t, ISO1_AIR_over_AUV, marker='x', color='black')
	ISO2 = plt.scatter(frc_t, ISO2_AIR_over_AUV, marker='x', color='red')
	ISO3 = plt.scatter(frc_t, ISO3_AIR_over_AUV, marker='x', color='blue')
	
	

	#plt.axvline(x=frc_PG, color='black', linestyle=':')
	plt.axvspan(frc_PG, frc_PGb, color='grey', alpha=0.5, lw=0)
	plt.axhline(y=0, color='black', linestyle='--')
	

	plt.legend( [ Anl1[0],  ISO1,  ISO2, ISO3], ("Analytic", r"$0.5\nu_0$", r"$\nu_0$",r"$2\nu_0$"), loc='upper right')

	plt.ylabel(r"$A_{\rm{IR}} / A$")
	plt.xlabel(r"$t_d / P$")

	

	plt.xlim(0.0,3.0)
	#plt.show()

	Savename = "plots/Iso_and_Dop/Analytics/Qv_varnu0_AIRoAUV_J%g_numin%g_numx%g_reclim2.png" %(JJt, Nnumn, Nnumx)
	Savename = Savename.replace('.', 'p')
	Savename = Savename.replace('ppng', '.png')
	plt.savefig(Savename)


