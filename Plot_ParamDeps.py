#TODO 3/18/16
#ADD GRIDLINES
#ADD LEGENDS
#
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
Thin  = True

Plot_I   = False
I_name = "Incs"
Plot_R   = False
R_name = "Rdust"
Plot_Om  = False
Om_name = "Ombins"
Plot_bet = False
bet_name = "betas"
Plot_J = True
J_name = "Js"
Plot_TT = True
TT_name = "Theta_Ts"


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
### DUST stuff
## for Qv
nne = 1.
nu0 = numicron*0.2
Rde = Rorb
Rrout = 60.0*Rde
pp = 1.1
thetTst = 1.*np.pi/4
JJt =4.*np.pi/8
aeff = 0.1*10**(-4) #(0.1 micrometer is an average ISM dust grain size)
md = 10**(-14)
n0 = 0.000001#6.*10**5*Msun/md * 1./(4./3.*ma.pi*(Rrout**3 - Rde**3))

##BINARY STUFF
Lav = L0
betst = 0.10
Inc = ma.acos(0.07/betst)#0.*np.pi/4.
Ombn = Omb
alph = 1.1
Dst = 1.4*10**9*pc2cm









### WISE BAND + Observational STUFF
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3

nuVbnd = c/(545*10**(-7))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-20)*(W1mn + W1mx)/2
FW2Rel = 1.7187*10**(-20)*(W2mn + W2mx)/2



##TABULATE T's and RHSs
print "Creating look up tables"
NT = 10000
RHS_table = np.zeros(NT)
T_table = np.linspace(0.1, 3000., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)


### PLOT POINTS
Nt=20
tt = np.linspace(0., 2.,       Nt)*2*np.pi/Ombn


## INTEGRATION LIMTS FOR ALL nu
numn = 0.00001*numicron
numx = 5*numicron


###########------------------------------------###########
###########		Plot Parameter Dependences     ###########
###########------------------------------------###########

if (Thin):
	###########---------------------------###########
	###########--------THIN SPHERE--------###########
	###########---------------------------###########
	####-------Inclination-------####
	Inc1 = 0.0
	#nrm = 1.27392461  #hack
	if (Plot_I):
		Inc1 = 0.0
		Inc2 = ma.pi/4.
		Inc3 = ma.pi/2.
		argI1 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		argI2 = [Lav, betst, Inc2, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		argI3 = [Lav, betst, Inc3, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
		argI1 = np.array(argI1)
		argI2 = np.array(argI2)
		argI3 = np.array(argI3)

		FsrcI1 = np.empty(Nt)
		FsrcI2 = np.empty(Nt)
		FsrcI3 = np.empty(Nt)
		FI1 = np.empty(Nt)
		FI2 = np.empty(Nt)
		FI3 = np.empty(Nt)

		
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FsrcI2 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc2, Ombn, alph)/FVbndRel)
		FsrcI3 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc3, Ombn, alph)/FVbndRel)
		#for i in range (0, Nt):
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argI1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argI2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argI3, RHS_table, T_table)/FW1Rel)

		nrm = FsrcI3 - FI3
		###PLOT###
		plt.figure()
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
		plt.xlim(0.0, 2.0)

		#plt.show()
		plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+I_name+"n0_%g.png" %n0)

	####-------END Inclination-------####





	####-------Rdust-------####
	Rin1 = Rorb
	if (Plot_R):
		Rin2 = Rorb/2.
		Rin3 = Rorb/4.
		argR1 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
		argR2 = [Lav, betst, Inc1, Ombn, alph, n0, Rin2, pp, thetTst, JJt, aeff, nu0, nne]
		argR3 = [Lav, betst, Inc1, Ombn, alph, n0, Rin3, pp, thetTst, JJt, aeff, nu0, nne]

		FsrcI1 = np.empty(Nt)
		FsrcI2 = np.empty(Nt)
		FsrcI3 = np.empty(Nt)
		FI1 = np.empty(Nt)
		FI2 = np.empty(Nt)
		FI3 = np.empty(Nt)

		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argR1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argR2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argR3, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F_{\rm{Bol}}$',   r'$R_d=R_{\rm{orb}}/c$',  r'$R_d=(R_{\rm{orb}}/c)/2$', r'$R_d=(R_{\rm{orb}}/c)/4$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(0.0, 2.0)

		#plt.show()
		plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+R_name+"n0_%g.png" %n0)


	####-------END Rdust-------####





	####-------Ombin-------####
	Om1 = Ombn
	if (Plot_Om):
		Om2 = Ombn*3.
		Om3 = Ombn*4.
		tt1 = np.linspace(0., 2.,       Nt)*2*np.pi/Om1
		tt2 = np.linspace(0., 2.,       Nt)*2*np.pi/Om2
		tt3 = np.linspace(0., 2.,       Nt)*2*np.pi/Om3
		argO1 = [Lav, betst, Inc1, Om1, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
		argO2 = [Lav, betst, Inc1, Om2, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
		argO3 = [Lav, betst, Inc1, Om3, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]

		FsrcI1 = np.empty(Nt)
		FsrcI2 = np.empty(Nt)
		FsrcI3 = np.empty(Nt)
		FI1 = np.empty(Nt)
		FI2 = np.empty(Nt)
		FI3 = np.empty(Nt)

		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt[i], Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, tt1, Dst, Rrout, argO1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, tt2, Dst, Rrout, argO2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, tt3, Dst, Rrout, argO3, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1) 
		###PLOT###
		plt.figure()
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)


		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F_{\rm{Bol}}$',   r'$\Omega=\Omega_{0}$',  r'$\Omega=3\Omega_{0}$', r'$\Omega=4\Omega_{0}$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(0.0, 2.0)

		#plt.show()
		plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+Om_name+"n0_%g.png" %n0)

	####-------END Ombin-------####




	####-------beta-------####
	bet1 = 0.1
	if (Plot_R):
			bet2 = 0.2
			bet3 = 0.3
			argb1 = [Lav, bet1, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
			argb2 = [Lav, bet2, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]
			argb3 = [Lav, bet3, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, JJt, aeff, nu0, nne]

			FsrcI1 = np.empty(Nt)
			FsrcI2 = np.empty(Nt)
			FsrcI3 = np.empty(Nt)
			FI1 = np.empty(Nt)
			FI2 = np.empty(Nt)
			FI3 = np.empty(Nt)

			#for i in range (0, Nt):
			FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, bet1, Inc1, Ombn, alph)/FVbndRel)
			FsrcI2 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, bet2, Inc1, Ombn, alph)/FVbndRel)
			FsrcI3 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, bet3, Inc1, Ombn, alph)/FVbndRel)
			FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argb1, RHS_table, T_table)/FW1Rel)
			FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argb2, RHS_table, T_table)/FW1Rel)
			FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argb3, RHS_table, T_table)/FW1Rel)

			nrm = np.mean(FsrcI1) - np.mean(FI1)
			###PLOT###
			plt.figure()

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
			plt.xlim(0.0, 2.0)

			#plt.show()
			plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+bet_name+"n0_%g.png" %n0)


	####-------END beta -------####



	####-------JJ-------####
	Rin1 = Rorb
	Inc1=0.0
	if (Plot_J):
		J1 = 0.
		J2 = ma.pi/2. - thetTst
		J3 = ma.pi/3.
		J4 = ma.pi/2. 
		argJ1 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J1, aeff, nu0, nne]
		argJ2 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J2, aeff, nu0, nne]
		argJ3 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J3, aeff, nu0, nne]
		argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J4, aeff, nu0, nne]

	

		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ2, RHS_table, T_table)/FW1Rel)
		#FI3 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ3, RHS_table, T_table)/FW1Rel)
		FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		#IR3=plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)
		IR4=plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='yellow', linewidth=2)

		plt.grid(b=True, which='both')
		#plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0], IR4[0] ], (r'$F_{\rm{Bol}}$', r'$J=0$',   r'$J=\pi/4$',  r'$J=\pi/3$', r'$J=\pi/2$'), loc='upper right')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR4[0] ], (r'$F_{\rm{Bol}}$', r'$J=0$',   r'$J=\pi/4$', r'$J=\pi/2$'), loc='upper right')


		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(0.0, 2.0)

		#plt.show()
		plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/"+J_name+"n0_%g.png" %n0)


	####-------END Rdust-------####


###########---------------------------###########
###########--------THICK TORUS--------###########
###########---------------------------###########
if (Thick):
	if (Plot_J):
		Inc1 = 0.0
		J1 = 0.
		J2 = ma.pi/2. - thetTst
		J3 = ma.pi/2.
		#J4 = ma.pi/2. + thetTst
		argJ1 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, J1, aeff, nu0, nne]
		argJ2 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, J2, aeff, nu0, nne]
		argJ3 = [Lav, betst, Inc1, Ombn, alph, n0, Rde, pp, thetTst, J3, aeff, nu0, nne]
		#argJ4 = [Lav, betst, Inc1, Ombn, alph, n0, Rin1, pp, thetTst, J4, aeff, nu0, nne]



		#for i in range (0, Nt):
		FsrcI1 = -2.5*np.log10(Fsrc(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
		FI1 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argJ1, RHS_table, T_table)/FW1Rel)
		FI2 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argJ2, RHS_table, T_table)/FW1Rel)
		FI3 = -2.5*np.log10(Fobs_Thick(numn, numx, tt, Dst, Rrout, argJ3, RHS_table, T_table)/FW1Rel)
		#FI4 = -2.5*np.log10(Fobs_Shell(numn, numx, tt, Dst, Rrout, argJ4, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3, color='brown', linewidth=2)
		#IR4=plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0] ], (r'$F_{\rm{Bol}}$', r'$J=0$',  r'$J=\pi/2 - \theta_T$', r'$J=\pi/2$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(0.0, 2.0)

		#plt.show()
		plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+J_name+"n0_%g.png" %n0)




	if (Plot_TT):
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

		nrm = np.mean(FsrcI1) - np.mean(FI1)
		###PLOT###
		plt.figure()
		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1, color='red', linewidth=2)
		s1=plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR2=plt.plot(tt/(2*np.pi/Ombn), FI2, color='orange', linewidth=2)

		IR3=plt.plot(tt/(2*np.pi/Ombn), FI3, color='brown', linewidth=2)
		#IR4=plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='brown', linewidth=2)

		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0] ], (r'$F_{\rm{Bol}}$', r'$\theta_T=0$',  r'$\theta_T=\pi/4$', r'$\theta_T=\pi/3$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(0.0, 2.0)

		#plt.show()
		plt.savefig("/Users/dorazio/Desktop/Current_Projects/MBHB_LightEchoes/python/Plot_ParamDep/Thick"+TT_name+"n0_%g.png" %n0)

