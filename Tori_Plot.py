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
FISO = True
FDOP = True

ShTor_OptThin = True
ShTor_OptThick = True
Thick_Tor = False





##ISO VARS:
Amp = 0.22105  ##For alpha=0.0, beta=0.0677
t0 = ma.pi/2. #0.0




#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
Omb = 1./(1*yr2sec)
L0 = 6.78*10**46
MPGmx = 10**9.4*Msun
Ryr = c*yr2sec
RdPG = ma.sqrt(0.1)*2.8 *pc2cm
OmPG = Omb*2.*ma.pi/4.1
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
thetTst = 1.*ma.pi/4.
JJt = 0.0#3*ma.pi/4.
aeff = 0.16*10**(-4) #(0.1 micrometer is an average ISM dust grain size - choose 0.16 to make nu0~1um)
md = 10**(-14)
n10 = 1.0/(ma.pi*Rde*aeff*aeff) * (pp-1.)
nfac = 1.2
n0 = nfac*n10 ##6.*10**5*Msun/md * 1./(4./3.*ma.pi*(Rrout**3 - Rde**3))


Lav = L0
betst = 0.06776  ## gives 0.14 mag amplitdue at I=0 
Inc = 0.0#ma.pi/2#ma.acos(0.07/betst)#0.*np.pi/4.
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
Nt=40
tt = np.linspace(0., 2.,       Nt)*2*ma.pi/Ombn


## INTEGRATION LIMTS FOR ALL nu
Nnumn = 0.00001
Nnumx = 100.0
numn = Nnumn*numicron
numx = Nnumx*numicron

# TD = np.zeros(Nt)
# arg1 = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
# for i in range(Nt):
# 	TD[i] = TDust_Iso(tt[i], Rde, ma.pi/2, 0.0, arg1, RHS_table, T_table)

# plt.figure()
# plt.plot(tt, TD)
# plt.show()
if (FISO):
###########------------------------------------###########
###########					 FISO	           ###########
###########------------------------------------###########
	if (ShTor_OptThin):
		print "Plottin F_ISO Shell Torus, Opt Thin"
		###########---------------------------###########
		###########-------------J-----------###########
		###########---------------------------###########

		print "Vary J"
		JJ1 = 0.0
		JJ2 = (thetTst - ma.pi/2)/2.
		JJ3 = (thetTst - ma.pi/2)
		JJ4 = ma.pi/2

		#Omfac = 10.
		#Ombn = Ombn*Omfac
		#tt = np.linspace(0., 2.,       Nt)*2*ma.pi/Ombn

		arg1 = [Lav, Amp, Ombn, t0, n0, 1.33*Rde, pp, thetTst, JJ1, aeff, nu0, nne]
		arg2 = [Lav, Amp, Ombn, t0, n0, 1.33*Rde, pp, thetTst, JJ2, aeff, nu0, nne]
		arg3 = [Lav, Amp, Ombn, t0, n0, 1.33*Rde, pp, thetTst, JJ3, aeff, nu0, nne]
		arg4 = [Lav, Amp, Ombn, t0, n0, 1.33*Rde, pp, thetTst, JJ4, aeff, nu0, nne]

		FI1 = np.zeros(Nt)
		FI2 = np.zeros(Nt)
		FI3 = np.zeros(Nt)
		FI4 = np.zeros(Nt)

		FsrcI1 = -2.5*np.log10(Fsrc_Iso(tt, Dst, Lav, Amp, Ombn, t0)/FVbndRel)
							#F_Ring_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
		for i in range(Nt):
			FI1[i] = -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
			FI2[i] = -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
			FI3[i] = -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)
			FI4[i] = -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt(numn, numx, tt[i], Dst, arg4, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1)

		###PLOT###
		plt.figure()
		#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
		NRd = Rde/3.08e18
		plt.title(r"Shell Torus A"+"\n"+r"$F^{\rm{src}}_{\rm{iso}}$, $\theta_T = \pi/4$")#" $\Omega = %i \times 2\pi c / R_d$" %(Omfac))

		s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)

		IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)

		IR4 = plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='yellow', linewidth=2)


		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0], IR4[0]  ], (r'$F^{\rm{src}}_{\rm{iso}}$', r'$J=0$',   r'$J=(\theta_T-\pi/2)/2$',   r'$J=(\theta_T-\pi/2)$',  r'$J=\pi/2$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
		#plt.ylim(plt.ylim(11.7, 12.7)[::-1])

		Savename = "plots/Iso_and_Dop/ShTor_Thin/FISO_ShTor_Thin_nrm%g_"%nrm+"_Rin%g_Inc%g_thetaT%g_VaryJ_numin%g_numx%g.png" %(Rde, Inc, thetTst, Nnumn, Nnumx)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)

		####-------END J-------####




	if (ShTor_OptThick):
		print "Plottin F_ISO Shell Torus, Opt Thick"
		###########---------------------------###########
		###########-------------J-----------###########
		###########---------------------------###########

		print "Vary J"
		JJ1 = 0.0
		JJ2 = (thetTst - ma.pi/2)/2.
		JJ3 = (thetTst - ma.pi/2)
		JJ4 = ma.pi/2

		#Omfac = 10.
		#Ombn = Ombn*Omfac
		#tt = np.linspace(0., 2.,       Nt)*2*ma.pi/Ombn

		arg1 = [Lav, Amp, Ombn, t0, n0, 1.33*Rde, pp, thetTst, JJ1, aeff, nu0, nne]
		arg2 = [Lav, Amp, Ombn, t0, n0, 1.33*Rde, pp, thetTst, JJ2, aeff, nu0, nne]
		arg3 = [Lav, Amp, Ombn, t0, n0, 1.33*Rde, pp, thetTst, JJ3, aeff, nu0, nne]
		arg4 = [Lav, Amp, Ombn, t0, n0, 1.33*Rde, pp, thetTst, JJ4, aeff, nu0, nne]

		FI1 = np.zeros(Nt)
		FI2 = np.zeros(Nt)
		FI3 = np.zeros(Nt)
		FI4 = np.zeros(Nt)

		FsrcI1 = -2.5*np.log10(Fsrc_Iso(tt, Dst, Lav, Amp, Ombn, t0)/FVbndRel)
							#F_Ring_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
		for i in range(Nt):
			FI1[i] = -2.5*np.log10(F_ShTorOptThick_Iso_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
			FI2[i] = -2.5*np.log10(F_ShTorOptThick_Iso_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
			FI3[i] = -2.5*np.log10(F_ShTorOptThick_Iso_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)
			FI4[i] = -2.5*np.log10(F_ShTorOptThick_Iso_QuadInt(numn, numx, tt[i], Dst, arg4, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1)

		###PLOT###
		plt.figure()
		#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
		NRd = Rde/3.08e18
		plt.title(r"Shell Torus B"+"\n"+r"$F^{\rm{src}}_{\rm{iso}}$, $\theta_T = \pi/4$")#" $\Omega = %i \times 2\pi c / R_d$" %(Omfac))

		s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)

		IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)

		IR4 = plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='yellow', linewidth=2)


		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0], IR4[0]  ], (r'$F^{\rm{src}}_{\rm{iso}}$', r'$J=0$',   r'$J=(\theta_T-\pi/2)/2$',   r'$J=(\theta_T-\pi/2)$',  r'$J=\pi/2$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
		#plt.ylim(plt.ylim(11.7, 12.7)[::-1])

		Savename = "plots/Iso_and_Dop/ShTor_Thick/FISO_ShTor_Thick_nrm%g_"%nrm+"_Rin%g_Inc%g_thetaT%g_VaryJ_numin%g_numx%g.png" %(Rde, Inc, thetTst, Nnumn, Nnumx)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)

		####-------END J-------####







if (FDOP):
###########------------------------------------###########
###########				DOPPLER		           ###########
###########------------------------------------###########
	print "Plottin F_Dop Shell Torus, Opt Thin"
	###########---------------------------###########
	###########-------------J-----------###########
	###########---------------------------###########
	if (ShTor_OptThin):
		print "Opt Thin - Vary J"
		JJ1 = 0.0
		JJ2 = (thetTst - ma.pi/2)/2.
		JJ3 = (thetTst - ma.pi/2)
		JJ4 = ma.pi/2

		#Omfac = 1.
		#Ombn = Ombn*Omfac
		#tt = np.linspace(0., 2.,       Nt)*2*ma.pi/Ombn

		arg1 = [Lav, betst, Inc, Ombn, alph, n0, 1.33*Rde, pp, thetTst, JJ1, aeff, nu0, nne]
		arg2 = [Lav, betst, Inc, Ombn, alph, n0, 1.33*Rde, pp, thetTst, JJ2, aeff, nu0, nne]
		arg3 = [Lav, betst, Inc, Ombn, alph, n0, 1.33*Rde, pp, thetTst, JJ3, aeff, nu0, nne]
		arg4 = [Lav, betst, Inc, Ombn, alph, n0, 1.33*Rde, pp, thetTst, JJ4, aeff, nu0, nne]

		FI1 = np.zeros(Nt)
		FI2 = np.zeros(Nt)
		FI3 = np.zeros(Nt)
		FI4 = np.zeros(Nt)

		FsrcI1 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)/FVbndRel)
							#F_Ring_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
		for i in range(Nt):
			FI1[i] = -2.5*np.log10(F_ShTorOptThin_Dop_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
			FI2[i] = -2.5*np.log10(F_ShTorOptThin_Dop_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
			FI3[i] = -2.5*np.log10(F_ShTorOptThin_Dop_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)
			FI4[i] = -2.5*np.log10(F_ShTorOptThin_Dop_QuadInt(numn, numx, tt[i], Dst, arg4, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1)

		###PLOT###
		plt.figure()
		#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
		NRd = Rde/3.08e18
		plt.title(r"Shell Torus A"+"\n"+r"$F^{\rm{src}}_{\rm{Dop}}$, $\theta_T = \pi/4$")#" $\Omega = %i \times 2\pi c / R_d$" %(Omfac))

		s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)

		IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)

		IR4 = plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='yellow', linewidth=2)


		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0], IR4[0]  ], (r'$F^{\rm{src}}_{\rm{Dop}}$', r'$J=0$',   r'$J=(\theta_T-\pi/2)/2$',   r'$J=(\theta_T-\pi/2)$',  r'$J=\pi/2$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
		#plt.ylim(plt.ylim(11.7, 12.7)[::-1])

		Savename = "plots/Iso_and_Dop/ShTor_Thin/FDop_ShTor_Thin_nrm%g_"%nrm+"_Rin%g_Inc%g_thetaT%g_VaryJ_numin%g_numx%g.png" %(Rde, Inc, thetTst, Nnumn, Nnumx)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)

		####-------END J-------####




	if (ShTor_OptThick):
		print "Plottin F_Dop Shell Torus, Opt Thick"
		###########---------------------------###########
		###########-------------J-----------###########
		###########---------------------------###########
		print "Opt Thick - Vary J"
		JJ1 = 0.0
		JJ2 = (thetTst - ma.pi/2)/2.
		JJ3 = (thetTst - ma.pi/2)
		JJ4 = ma.pi/2

		#Omfac = 10.
		#Ombn = Ombn*Omfac
		#tt = np.linspace(0., 2.,       Nt)*2*ma.pi/Ombn

		arg1 = [Lav, betst, Inc, Ombn, alph, n0, 1.33*Rde, pp, thetTst, JJ1, aeff, nu0, nne]
		arg2 = [Lav, betst, Inc, Ombn, alph, n0, 1.33*Rde, pp, thetTst, JJ2, aeff, nu0, nne]
		arg3 = [Lav, betst, Inc, Ombn, alph, n0, 1.33*Rde, pp, thetTst, JJ3, aeff, nu0, nne]
		arg4 = [Lav, betst, Inc, Ombn, alph, n0, 1.33*Rde, pp, thetTst, JJ4, aeff, nu0, nne]

		FI1 = np.zeros(Nt)
		FI2 = np.zeros(Nt)
		FI3 = np.zeros(Nt)
		FI4 = np.zeros(Nt)

		FsrcI1 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)/FVbndRel)
							#F_Ring_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
		for i in range(Nt):
			FI1[i] = -2.5*np.log10(F_ShTorOptThick_Dop_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
			FI2[i] = -2.5*np.log10(F_ShTorOptThick_Dop_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
			FI3[i] = -2.5*np.log10(F_ShTorOptThick_Dop_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)
			FI4[i] = -2.5*np.log10(F_ShTorOptThick_Dop_QuadInt(numn, numx, tt[i], Dst, arg4, RHS_table, T_table)/FW1Rel)

		nrm = np.mean(FsrcI1) - np.mean(FI1)

		###PLOT###
		plt.figure()
		#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
		NRd = Rde/3.08e18
		plt.title(r"Shell Torus B"+"\n"+r"$F^{\rm{src}}_{\rm{Dop}}$, $\theta_T = \pi/4$")#" $\Omega = %i \times 2\pi c / R_d$" %(Omfac))

		s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=2)

		IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='red', linewidth=2)

		IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='orange', linewidth=2)

		IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='brown', linewidth=2)

		IR4 = plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='yellow', linewidth=2)


		plt.grid(b=True, which='both')
		plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0], IR4[0]  ], (r'$F^{\rm{src}}_{\rm{Dop}}$', r'$J=0$',   r'$J=(\theta_T-\pi/2)/2$',   r'$J=(\theta_T-\pi/2)$',  r'$J=\pi/2$'), loc='upper right')

		plt.xlabel(r"$N_{\rm{orb}}$")
		plt.ylabel("mag")
		plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
		#plt.ylim(plt.ylim(11.7, 12.7)[::-1])

		Savename = "plots/Iso_and_Dop/ShTor_Thick/FDop_ShTor_Thick_nrm%g_"%nrm+"_Rin%g_Inc%g_thetaT%g_VaryJ_numin%g_numx%g.png" %(Rde, Inc, thetTst, Nnumn, Nnumx)
		Savename = Savename.replace('.', 'p')
		Savename = Savename.replace('ppng', '.png')
		plt.savefig(Savename)

		####-------END J-------####
