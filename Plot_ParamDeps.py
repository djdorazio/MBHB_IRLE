import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] = 'sans-serif'
#matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 16})

import IR_LightEchoes_NewMeth as IRLE
from IR_LightEchoes_NewMeth import *






###OPTIONS
Thick = True
Thin  = False
Plot_v_R = False


Plot_I   = False
I_name = "Incs"
Plot_R   = False
R_name = "Rdust"
Plot_Om  = False
Om_name = "Ombins"
Plot_bet = False
bet_name = "betas"
Plot_J = False
J_name = "Js"
Plot_TT = False
TT_name = "Theta_Ts"
Plot_pp = False
pp_name = "_ps"
Plot_Ro = False
Ro_name = "_Routs"
Plot_n0 = True


#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
Omb = 1./(1*yr2sec)
L0 = 6.78*10**46
MPGmx = 10**9.4*Msun
Ryr = c*yr2sec
RdPG = ma.sqrt(0.1)*2.8 *pc2cm
OmPG = Omb*2.*np.pi/4.1
#alphnu = 1.1

Rorb = c*2.*ma.pi/Omb
Ompc = 2.*ma.pi*c/pc2cm/2.


## TEST VALUES
### DUST stuff
## for Qv
nne = 1.
nu0 = numicron
Rde = RdPG
Rrout = 10.0*Rde
pp = 2.0
thetTst = 1.*np.pi/4.
JJt =np.pi/2.
aeff = 0.16*10**(-4) #(0.1 micrometer is an average ISM dust grain size - choose 0.16 to make nu0~1um)
md = 10**(-14)
n10 = 0.1/(ma.pi*Rde*aeff*aeff) * (pp-1.)
nfac = 10.0
n0 = nfac*n10 ##6.*10**5*Msun/md * 1./(4./3.*ma.pi*(Rrout**3 - Rde**3))



##BINARY STUFF
Lav = L0
betst = 0.10
Inc = 0.0#ma.acos(0.07/betst)#0.*np.pi/4.
Ombn = 2.*ma.pi/(Rde/c)  ##(2pi/P)
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



##TABULATE T's and RHSs
print "Creating look up tables"
NT = 10000
RHS_table = np.zeros(NT)
T_table = np.linspace(1., 2000., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)


### PLOT POINTS
Nt=20
tt = np.linspace(0., 2.,       Nt)*2*np.pi/Ombn


## INTEGRATION LIMTS FOR ALL nu
numn = 0.001*numicron
numx = 10*numicron


###########------------------------------------###########
###########		Plot Parameter Dependences     ###########
###########------------------------------------###########

if (Thin):
	print "Plottin Shells"
	###########---------------------------###########
	###########--------THIN SPHERE--------###########
	###########---------------------------###########
	####-------Inclination-------####
	Inc1 = 0.0
	if (Plot_I):
		print "Vary Inc Shells"
		Inc1 = 0.0
		Inc2 = ma.pi/4.
		Inc3 = ma.pi/2.
		argI1 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		argI2 = [Lav, betst, Inc2, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		argI3 = [Lav, betst, Inc3, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		# argI1 = np.array(argI1)
		# argI2 = np.array(argI2)
		# argI3 = np.array(argI3)

		# FsrcI1 = np.empty(Nt)
		# FsrcI2 = np.empty(Nt)
		# FsrcI3 = np.empty(Nt)
		# FI1 = np.empty(Nt)
		# FI2 = np.empty(Nt)
		# FI3 = np.empty(Nt)

		
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FsrcI2 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc2, Ombn, alph)/FVbndRel)
		FsrcI3 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc3, Ombn, alph)/FVbndRel)
		
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, Rde, tt, Dst, Rrout, argI1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, Rde, tt, Dst, Rrout, argI2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, Rde, tt, Dst, Rrout, argI3, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI3) - np.mean(FI3)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)
		s2 = plt.plot(tt/(2*np.pi/Ombn), FsrcI2, linestyle = '--', color='purple', linewidth=2)

		IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)
		s3 = plt.plot(tt/(2*np.pi/Ombn), FsrcI3, linestyle = '--', color='green', linewidth=2)

		
		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], s2[0], IR2[0], s3[0], IR3[0]  ], (r'$i=0$','',   r'$i=\pi/4$','',   r'$i=\pi/2$', ''), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))


		Savename = "plots/Shell_nrm%g_"%nrm+I_name+"_Rin%g_J%g_thetT%g_Rout%g_p%g_n0%g.png" %(Rde, JJt, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)

	####-------END Inclination-------####





	####-------Rdust-------####
	if (Plot_R):
		print "Vary Rin Shells"
		Rin1 = Rde
		Rin2 = Rde*2.
		Rin3 = Rde*3.
		argR1 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
		argR2 = [Lav, betst, Inc1, Ombn, alph, n0, Rin2, pp, thetTst, JJt, aeff, nu0, nne]
		argR3 = [Lav, betst, Inc1, Ombn, alph, n0, Rin3, pp, thetTst, JJt, aeff, nu0, nne]


		
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1    = -2.5*np.log10(Fobs_Shell(numn, numx, Rin1, tt, Dst, Rrout, argR1, RHS_table, T_table)/FW1Rel)
		FI2    = -2.5*np.log10(Fobs_Shell(numn, numx, Rin2, tt, Dst, Rrout, argR2, RHS_table, T_table)/FW1Rel)
		FI3    = -2.5*np.log10(Fobs_Shell(numn, numx, Rin3, tt, Dst, Rrout, argR3, RHS_table, T_table)/FW1Rel)

		nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.title(r"$n_0 = %g n_T$" %nfac)
		plt.figure()
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F_{\rm{Bol}}$',   r'$R_d=R_{\rm{orb}}/c$',  r'$R_d=2(R_{\rm{orb}}/c)$', r'$R_d=3(R_{\rm{orb}}/c)$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
		#plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+R_name+"n0_%g.png" %n0)
		#plt.savefig("plots/"+R_name+"n0_%g.png" %n0)
		Savename = "plots/Shell_nrm%g_"%nrm+R_name+"_Om%g_J%g_thetT%g_Rout%g_p%g_n0%g.png" %(Ombn, JJt, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)

	####-------END Rdust-------####





	####-------Ombin-------####
	
	if (Plot_Om):
		print "Vary Om Shells"
		Om1 = 0.5*Ombn
		Om2 = Ombn*1.
		Om3 = Ombn*1.5
		tt1 = np.linspace(0., 2.,       Nt)*2*np.pi/Om1
		tt2 = np.linspace(0., 2.,       Nt)*2*np.pi/Om2
		tt3 = np.linspace(0., 2.,       Nt)*2*np.pi/Om3
		argO1 = [Lav, betst, Inc1, Om1, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
		argO2 = [Lav, betst, Inc1, Om2, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
		argO3 = [Lav, betst, Inc1, Om3, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]

		# FsrcI1 = np.empty(Nt)
		# FsrcI2 = np.empty(Nt)
		# FsrcI3 = np.empty(Nt)
		# FI1 = np.empty(Nt)
		# FI2 = np.empty(Nt)
		# FI3 = np.empty(Nt)

		rem = Rde
		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, rem, tt1, Dst, Rrout, argO1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, rem, tt2, Dst, Rrout, argO2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, rem, tt3, Dst, Rrout, argO3, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1) 
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt1/(2*np.pi/Om1), FI1+nrm, color='red', linewidth=2)
		s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2 = plt.plot(tt2/(2*np.pi/Om2), FI2+nrm, color='orange', linewidth=2)

		IR3 = plt.plot(tt3/(2*np.pi/Om3), FI3+nrm, color='brown', linewidth=2)


		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F_{\rm{Bol}}$',   r'$\Omega=0.5\Omega_{0}$',  r'$\Omega=\Omega_{0}$', r'$\Omega=1.5\Omega_{0}$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
		#plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+Om_name+"n0_%g.png" %n0)
		#plt.savefig("plots/"+Om_name+"n0_%g.png" %n0)
		Savename = "plots/Shell_nrm%g_"%nrm+Om_name+"_Rin%g_J%g_thetT%g_Rout%g_p%g_n0%g.png" %(Rde, JJt, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
	####-------END Ombin-------####




	####-------beta-------####
	bet1 = 0.1
	if (Plot_bet):
		print "Vary beta Shells"
		bet2 = 0.2
		bet3 = 0.3
		argb1 = [Lav, bet1, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
		argb2 = [Lav, bet2, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
		argb3 = [Lav, bet3, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]

		# FsrcI1 = np.empty(Nt)
		# FsrcI2 = np.empty(Nt)
		# FsrcI3 = np.empty(Nt)
		# FI1 = np.empty(Nt)
		# FI2 = np.empty(Nt)
		# FI3 = np.empty(Nt)

		rem = Rde
		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, bet1, Inc1, Ombn, alph)/FVbndRel)
		FsrcI2 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, bet2, Inc1, Ombn, alph)/FVbndRel)
		FsrcI3 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, bet3, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, rem, tt, Dst, Rrout, argb1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, rem, tt, Dst, Rrout, argb2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, rem, tt, Dst, Rrout, argb3, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)
		s2 = plt.plot(tt/(2*np.pi/Ombn), FsrcI2, linestyle = '--', color='purple', linewidth=2)

		IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)
		s3 = plt.plot(tt/(2*np.pi/Ombn), FsrcI3, linestyle = '--', color='green', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], s2[0], IR2[0], s3[0], IR3[0]  ], (r'$\beta=0.1$','',   r'$\beta=0.2$','',   r'$\beta=0.3$', ''), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
		#plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+bet_name+"n0_%g.png" %n0)
		#plt.savefig("plots/"+bet_name+"n0_%g.png" %n0)
		Savename = "plots/Shell_nrm%g_"%nrm+bet_name+"_Rin%g_J%g_thetT%g_Rout%g_p%g_n0%g.png" %(Rde, JJt, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)

	####-------END beta -------####



	####-------JJ-------####
	Rin1 = Rde
	Inc1=0.0
	if (Plot_J):
		print "Vary J Shells"
		J1 = 0.
		J2 = ma.pi/2. - thetTst
		J3 = ma.pi/3.
		J4 = ma.pi/2. 
		argJ1 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J1, aeff, nu0, nne]
		argJ2 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J2, aeff, nu0, nne]
		argJ3 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J3, aeff, nu0, nne]
		argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J4, aeff, nu0, nne]

	
		rem = Rde
		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, rem, tt, Dst, Rrout, argJ1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, rem, tt, Dst, Rrout, argJ2, RHS_table, T_table)/FW1Rel)
		#FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ3, RHS_table, T_table)/FW1Rel)
		FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, rem, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		#IR3=plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)
		IR4=plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		#plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0], IR4[0] ], (r'$F_{\rm{Bol}}$', r'$J=0$',   r'$J=\pi/4$',  r'$J=\pi/3$', r'$J=\pi/2$'), loc='upper right')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR4[0] ], (r'$F_{\rm{Bol}}$', r'$J=0$',   r'$J=\pi/4$', r'$J=\pi/2$'), loc='upper right')


		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
		#plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+J_name+"n0_%g.png" %n0)
		#plt.savefig("plots/"+J_name+"n0_%g.png" %n0)
		Savename = "plots/Shell_nrm%g_"%nrm+J_name+"_Rin%g_thetT%g_Rout%g_p%g_n0%g.png" %(Rde, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
	####-------END JJ -------####



		####-------TT-------####
	if (Plot_TT):
		print "Vary theta_T Shells"
		Inc1 = 0.0
		JJ = ma.pi/2.
		TT1 = 0.
		TT2 = ma.pi/4
		TT3 = ma.pi/3.
		#J4 = ma.pi/2. + thetTst
		argJ1 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, TT1, JJ, aeff, nu0, nne]
		argJ2 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, TT2, JJ, aeff, nu0, nne]
		argJ3 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, TT3, JJ, aeff, nu0, nne]
		#argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J4, aeff, nu0, nne]



		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, Rde, tt, Dst, Rrout, argJ1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, Rde, tt, Dst, Rrout, argJ2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, Rde, tt, Dst, Rrout, argJ3, RHS_table, T_table)/FW1Rel)
		#FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3, color='brown', linewidth=2)
		#IR4=plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0] ], (r'$F_{\rm{Bol}}$', r'$\theta_T=0$',  r'$\theta_T=\pi/4$', r'$\theta_T=\pi/3$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
	#	plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+TT_name+"n0_%g.png" %n0)
		Savename = "plots/Shell_nrm%g_"%nrm+TT_name+"_J%g_Rout%g_p%g_n0%g.png" %(JJt, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
		####-------END TT-------####

		####-------Rout Shell-------####
	if (Plot_Ro):
		print "Vary Rout Shells"
		Inc1 = 0.0
		JJ = ma.pi/2.
		Ro1 = 10.*Rde
		Ro2 = 20.*Rde
		Ro3 = 100.*Rde
		#J4 = ma.pi/2. + thetTst
		argRo1 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, TT1, JJ, aeff, nu0, nne]
		argRo2 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, TT2, JJ, aeff, nu0, nne]
		argRo3 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, TT3, JJ, aeff, nu0, nne]
		#argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J4, aeff, nu0, nne]



		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, Rde, tt, Dst, Ro1, argJ1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, Rde, tt, Dst, Ro2, argJ2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, Rde, tt, Dst, Ro3, argJ3, RHS_table, T_table)/FW1Rel)
		#FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3, color='brown', linewidth=2)
		#IR4=plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0] ], (r'$F_{\rm{Bol}}$', r'$R_{\rm{out}}=10 R_{\rm{d}}$',  r'$R_{\rm{out}}=20 R_{\rm{d}}$', r'$R_{\rm{out}}=100 R_{\rm{d}}$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
	#	plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+TT_name+"n0_%g.png" %n0)
		Savename = "plots/Shell_nrm%g_"%nrm+Ro_name+"_J%g_thetT%g_p%g_n0%g.png" %(JJt, thetTst, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
		####-------END Rout Shell-------####

###########---------------------------###########
###########--------THICK TORUS--------###########
###########---------------------------###########
if (Thick):
	print "Plottin Thick"
	####-------THICK INC-------####
	if (Plot_I):
		print "Vary Inc Thick"
		Inc1 = 0.0
		Inc2 = ma.pi/3
		Inc3 =  ma.pi/2.

		JJt = ma.pi/2.

		argI1 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		argI2 = [Lav, betst, Inc2, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		argI3 = [Lav, betst, Inc3, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]

		FsrcI1 = np.empty(Nt)
		FsrcI2 = np.empty(Nt)
		FsrcI3 = np.empty(Nt)
		FI1 = np.empty(Nt)
		FI2 = np.empty(Nt)
		FI3 = np.empty(Nt)

		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FsrcI2 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc2, Ombn, alph)/FVbndRel)
		FsrcI3 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc3, Ombn, alph)/FVbndRel)

		FI1    = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argI1, RHS_table, T_table)/FW1Rel)
		FI2    = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argI2, RHS_table, T_table)/FW1Rel)
		FI3    = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argI3, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI3) - np.mean(FI3)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJ))
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)
		s2=plt.plot(tt/(2*np.pi/Ombn), FsrcI2, linestyle = '--', color='purple', linewidth=2)


		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)
		s3=plt.plot(tt/(2*np.pi/Ombn), FsrcI3, linestyle = '--', color='green', linewidth=2)


		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], s2[0], IR2[0], s3[0], IR3[0]  ], (r'$I=0$','',   r'$I=\pi/4$','',   r'$I=\pi/2$', ''), loc='upper right')
		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
		#plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+R_name+"n0_%g.png" %n0)
		Savename = "plots/Thick_nrm%g_"%nrm+I_name+"_J%g_thetT%g_Rout%g_p%g_n0%g.png" %(JJt, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)	
	####-------END THICK INC-------####

	####-------Rdust-------####
	if (Plot_R):
		print "Vary Rin Thick"
		Inc1 = 0.0
		Rin1 = Rde
		Om1 = 2.* ma.pi*c/Rde
		Rin2 = Rde*2.
		Om2 = 2.* ma.pi*c/(Rin2)
		Rin3 = Rde*3.
		Om3 = 2.* ma.pi*c/(Rin3)

		tt1 = np.linspace(0., 2.,       Nt)*2*np.pi/Om1
		tt2 = np.linspace(0., 2.,       Nt)*2*np.pi/Om2
		tt3 = np.linspace(0., 2.,       Nt)*2*np.pi/Om3

		argR1 = [Lav, betst, Inc1, Om1, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
		argR2 = [Lav, betst, Inc1, Om1, alph, n0, Rin2, pp, thetTst, JJt, aeff, nu0, nne]
		argR3 = [Lav, betst, Inc1, Om1, alph, n0, Rin3, pp, thetTst, JJt, aeff, nu0, nne]

		FsrcI1 = np.empty(Nt)
		FsrcI2 = np.empty(Nt)
		FsrcI3 = np.empty(Nt)
		FI1 = np.empty(Nt)
		FI2 = np.empty(Nt)
		FI3 = np.empty(Nt)

		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Om1, alph)/FVbndRel)
		FI1    = -2.5*np.log10(Fobs_Thick(numn, numx, tt1, Dst, Rrout, argR1, RHS_table, T_table)/FW1Rel)
		FI2    = -2.5*np.log10(Fobs_Thick(numn, numx, tt1, Dst, Rrout, argR2, RHS_table, T_table)/FW1Rel)
		FI3    = -2.5*np.log10(Fobs_Thick(numn, numx, tt1, Dst, Rrout, argR3, RHS_table, T_table)/FW1Rel)

		nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt/(2*np.pi/Om1), FI1+nrm, color='red', linewidth=2)
		s1=plt.plot(tt1/(2*np.pi/Om1), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt1/(2*np.pi/Om1), FI2+nrm, color='orange', linewidth=2)

		IR3=plt.plot(tt1/(2*np.pi/Om1), FI3+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F_{\rm{Bol}}$',   r'$R_d=R_{\rm{orb}}/c$',  r'$R_d=2(R_{\rm{orb}}/c)$', r'$R_d=3(R_{\rm{orb}}/c)$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
		#plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+R_name+"n0_%g.png" %n0)
		Savename = "plots/Thick_nrm%g_"%nrm+R_name+"_Om1_J%g_thetT%g_Rout%g_p%g_n0%g.png" %(JJt, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
	####-------END Rdust-------####

	####-------Ombin-------####
	if (Plot_Om):
		print "Vary Om Thick"
		Om1 = 0.5*Ombn
		Inc1 = 0.0
		Om2 = 1.0*Ombn
		Om3 = 1.5*Ombn
		tt1 = np.linspace(0., 2.,       Nt)*2*np.pi/Om1
		tt2 = np.linspace(0., 2.,       Nt)*2*np.pi/Om2
		tt3 = np.linspace(0., 2.,       Nt)*2*np.pi/Om3
		argO1 = [Lav, betst, Inc1, Om1, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		argO2 = [Lav, betst, Inc1, Om2, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		argO3 = [Lav, betst, Inc1, Om3, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]

		FsrcI1 = np.empty(Nt)
		FsrcI2 = np.empty(Nt)
		FsrcI3 = np.empty(Nt)
		FI1 = np.empty(Nt)
		FI2 = np.empty(Nt)
		FI3 = np.empty(Nt)

		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Thick(numn, numx, tt1, Dst, Rrout, argO1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Thick(numn, numx, tt2, Dst, Rrout, argO2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Thick(numn, numx, tt3, Dst, Rrout, argO3, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1) 
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt1/(2*np.pi/Om1), FI1+nrm, color='red', linewidth=2)
		s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2 = plt.plot(tt2/(2*np.pi/Om2), FI2+nrm, color='orange', linewidth=2)

		IR3 = plt.plot(tt3/(2*np.pi/Om3), FI3+nrm, color='brown', linewidth=2)


		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F_{\rm{Bol}}$',   r'$\Omega=0.5\Omega_{0}$',  r'$\Omega=\Omega_{0}$', r'$\Omega=1.5\Omega_{0}$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
		#plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+Om_name+"n0_%g.png" %n0)
		Savename = "plots/Thick_nrm%g_"%nrm+Om_name+"_J%g_thetT%g_Rout%g_p%g_n0%g.png" %(JJt, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
	####-------END Ombin-------####

	####-------JJ-------####
	if (Plot_J):
		print "Vary J Thick"
		Inc1 = 0.0
		J1 = 0.
		J2 = (ma.pi/2. - thetTst)
		J3 = ma.pi/2.
		#J4 = -(ma.pi/2. - thetTst)


		#J4 = ma.pi/2. + thetTst
		argJ1 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, J1, aeff, nu0, nne]
		argJ2 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, J2, aeff, nu0, nne]
		argJ3 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, J3, aeff, nu0, nne]
		#argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, J4, aeff, nu0, nne]



		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argJ1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argJ2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argJ3, RHS_table, T_table)/FW1Rel)
		#FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3, color='brown', linewidth=2)
		#IR4=plt.plot(tt/(2*np.pi/Ombn), FI4, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		#plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0], IR4[0] ], (r'$F_{\rm{Bol}}$', r'$J=0$',  r'$J=\pi/2 - \theta_T$', r'$J=\pi/2$', r'$J=-(\pi/2 - \theta_T)$'), loc='upper right')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0] ], (r'$F_{\rm{Bol}}$', r'$J=0$',  r'$J=\pi/2 - \theta_T$', r'$J=\pi/2$'), loc='upper right')


		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
	#	plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+J_name+"n0_%g.png" %n0)
		Savename = "plots/Thick_nrm%g_"%nrm+J_name+"_thetT%g_Rout%g_p%g_n0%g.png" %(thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
	####-------END JJ-------####

	####-------TT-------####
	if (Plot_TT):
		print "Vary theta_T Thick"
		Inc1 = 0.0
		JJ = ma.pi/2.
		TT1 = 0.
		TT2 = ma.pi/4
		TT3 = ma.pi/3.
		#J4 = ma.pi/2. + thetTst
		argJ1 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, TT1, JJ, aeff, nu0, nne]
		argJ2 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, TT2, JJ, aeff, nu0, nne]
		argJ3 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, TT3, JJ, aeff, nu0, nne]
		#argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J4, aeff, nu0, nne]



		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argJ1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argJ2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argJ3, RHS_table, T_table)/FW1Rel)
		#FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3, color='brown', linewidth=2)
		#IR4=plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0] ], (r'$F_{\rm{Bol}}$', r'$\theta_T=0$',  r'$\theta_T=\pi/4$', r'$\theta_T=\pi/3$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
	#	plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+TT_name+"n0_%g.png" %n0)
		Savename = "plots/Thick_nrm%g_"%nrm+TT_name+"_J%g_Rout%g_p%g_n0%g.png" %(JJt, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
		####-------END TT-------####



		####-------pp-------####
	if (Plot_pp):
		print "Vary p Thick"
		Inc1 = 0.0
		JJ = ma.pi/2.
		pp1 = 2.0
		pp2 = 3.0
		pp3 = 4.0
		#J4 = ma.pi/2. + thetTst
		argp1 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp1, thetTst, JJ, aeff, nu0, nne]
		argp2 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp2, thetTst, JJ, aeff, nu0, nne]
		argp3 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp3, thetTst, JJ, aeff, nu0, nne]
		#argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J4, aeff, nu0, nne]



		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argp1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argp2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argp3, RHS_table, T_table)/FW1Rel)
		#FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)
		#IR4=plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0] ], (r'$F_{\rm{Bol}}$', r'$p=2$',  r'$p=3$', r'$p=4$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
	#	plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+TT_name+"n0_%g.png" %n0)
		Savename = "plots/Thick_nrm%g_"%nrm+pp_name+"_J%g_thetT%g_Rout%g_n0%g.png" %(JJt, thetTst, Rrout, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
		####-------END TT-------####


		####------Rout-------####
	if (Plot_Ro):
		print "Vary Rout Thick"
		Inc1 = 0.0
		JJ = ma.pi/2.
		Ro1 = 10.*Rde
		Ro2 = 20.0*Rde
		Ro3 = 100.0*Rde
		#J4 = ma.pi/2. + thetTst
		arg = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, JJ, aeff, nu0, nne]
		arg = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, JJ, aeff, nu0, nne]
		arg = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, JJ, aeff, nu0, nne]
		#argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J4, aeff, nu0, nne]



		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Ro1, arg, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Ro2, arg, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Ro3, arg, RHS_table, T_table)/FW1Rel)
		#FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)
		#IR4=plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0] ], (r'$F_{\rm{Bol}}$', r'$R_{\rm{out}}=10 R_{\rm{d}}$',  r'$R_{\rm{out}}=20 R_{\rm{d}}$', r'$R_{\rm{out}}=100 R_{\rm{d}}$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

		#plt.show()
	#	plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+TT_name+"n0_%g.png" %n0)
		Savename = "plots/Thick_nrm%g_"%nrm+Ro_name+"_J%g_thetT%g_p%g_n0%g.png" %(JJt, thetTst, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
		####-------END pp-------####

	####------v n0-------####
	if (Plot_n0):
		Inc1 = 0.0
		JJ = ma.pi/2.
		Ro1 = 10.*Rde
		#J4 = ma.pi/2. + thetTst
		n0_arr = np.logspace(-1, 2.0, Nt)
		arg = []
		for i in range(Nt):
			arg.append([Lav, 0.0, Inc1, Ombn, alph, 10**(n0_arr[i]) *n10, Rde, pp, thetTst, JJ, aeff, nu0, nne])
		#argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J4, aeff, nu0, nne]



		#for i in range (0, Nt):
		#FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Thick_n0(numn, numx, 0.0, Dst, Ro1, arg, RHS_table, T_table)/FW1Rel)
		#FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		#nrm = np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		IR1 = plt.plot(n0_arr, FI1, color='red', linewidth=2)
		

		plt.grid(b=True, which='both')
		#plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0] ], (r'$F_{\rm{Bol}}$', r'$R_{\rm{out}}=10 R_{\rm{d}}$',  r'$R_{\rm{out}}=20 R_{\rm{d}}$', r'$R_{\rm{out}}=100 R_{\rm{d}}$'), loc='upper right')

		plt.xlabel(r"$n_0/n_T$")
		plt.ylabel(r"$\left<\rm{mag}\right>$")
		plt.xlim(0.1, 100.0)

		#plt.show()
	#	plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+TT_name+"n0_%g.png" %n0)
		plt.savefig("plots/Thick_varyn0_in units_of_n0_%g.png" %n0)

		Savename = "plots/Thick_varyn0_in units_of_n0_Rin%g_J%g_thetT%g_Rout%g_p%g_n0%g.png" %(Rde, JJt, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
		####-------END v n0-------####

	####------v R-------####
if (Plot_v_R):
		print "plot Fnu v nu for diff Routs"
		Inc1 = 0.0
		#JJ = ma.pi/2.
		Rout1 = Rde*5.
		Rout2 = Rde*20.
		Rout3 = Rde*50.
		Rout4 = Rde*100.
		#J4 = ma.pi/2. + thetTst
		nu_arr = np.linspace(0.001*numicron, 1.*numicron, 30)

		JJt = ma.pi/2#-(ma.pi/2. - thetTst) # so we can see all of back inner edge
		#n0 = 100.0*n0
		arg = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		
		FI1 = Fnu_Thick_mult(nu_arr, 0.0, Dst, Rout1, arg, RHS_table, T_table)
		FI2 = Fnu_Thick_mult(nu_arr, 0.0, Dst, Rout2, arg, RHS_table, T_table)
		FI3 = Fnu_Thick_mult(nu_arr, 0.0, Dst, Rout3, arg, RHS_table, T_table)
		FI4 = Fnu_Thick_mult(nu_arr, 0.0, Dst, Rout4, arg, RHS_table, T_table)
		
		

		#nrm = np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		plt.title(r"$n_0 = %g n_T$" %nfac)
		IR1 = plt.plot(nu_arr, FI1, color='orange', linewidth=2)
		IR2 = plt.plot(nu_arr, FI2, color='red', linewidth=2)
		IR3 = plt.plot(nu_arr, FI3, color='brown', linewidth=2)
		IR4 = plt.plot(nu_arr, FI4, color='black', linewidth=2)
		
		plt.axvline(W1mn, linestyle="--", color="black")
		plt.axvline(W1mx, linestyle="--", color="black")

		plt.axvline(W2mn, linestyle="--", color="blue")
		plt.axvline(W2mx, linestyle="--", color="blue")

		plt.grid(b=True, which='both')
		plt.legend( [ IR1[0], IR2[0], IR3[0], IR4[0] ], ( r'$R_{\rm{out}}=5 R_{\rm{d}}$',  r'$R_{\rm{out}}=20 R_{\rm{d}}$', r'$R_{\rm{out}}=50 R_{\rm{d}}$', r'$R_{\rm{out}}=100 R_{\rm{d}}$'), loc='upper right')
	#	plt.legend( [ IR1[0], IR4[0] ], (r'$R_{\rm{out}}=10 R_{\rm{d}}$', r'$R_{\rm{out}}=100 R_{\rm{d}}$'), loc='upper right')


		plt.xlabel(r"$\nu$ [Hz]")
		plt.ylabel(r"$F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$]")
		plt.xlim(0.01*numicron, 1.*numicron)

		#plt.show()
	#	plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+TT_name+"n0_%g.png" %n0)
		
		Savename = "plots/Thick_varynu_diffRout_J%g_thetT%g_Rout%g_p%g_n0%g.png" %(JJt, thetTst, Rrout, pp, n0)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)
		####-------END v n0-------####



