import matplotlib
#matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
#matplotlib.rcParams['font.family'] = 'sans-serif'
#matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 20})


from FluxFuncs_IRLE import *






###OPTIONS###OPTIONS
#which model
FISO = True
FDOP = False
Sphere = True
Ring = False

#Which variable
Plot_R   = True  # Fig 6 ##with Sphere True (Both Iso and Dop)
Plot_J   = False
Plot_I   = False  # Fig 7##with Sphere True and FDOP True
# Is there any alpha dependence on phase on doppler case???
Plot_alp = False  ## Diagnostic for sign of AIR/A and shows that light curves flip aroudn alpha=3


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
Amp = 0.15 #For alpha = 1, bet = 0.07 #0.22105  ##For alpha=0.0, beta=0.0677
t0 = ma.pi/2. #0.0




#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
Omb = 1./(1*yr2sec)
L0 = 6.78*10**46
MPGmx = 10**9.4*Msun
Ryr = c*yr2sec
RdPG = 1.0*pc2cm#ma.sqrt(0.1)*2.8 *pc2cm  ##changed in final version
OmPG = Omb*2.*ma.pi/4.1
#alphnu = 1.1

Rorb = c*2.*ma.pi/Omb
Ompc = 2.*ma.pi*c/pc2cm/2.


## TEST VALUES
### DUST stuff
## for Qv
nne = 1.6 ##changed in final version
nu0 = numicron
Rde = RdPG
Rrout = 1.0*Rde
pp = 2.0
thetTst = 0.0#1.*ma.pi/2.
JJt = 0.0#3*ma.pi/4.
aeff = (c/nu0)/(2.*ma.pi) #(0.1 micrometer is an average ISM dust grain size - choose 0.16 to make nu0~1um)
#md = 10**(-14)
n10 = 1.0/(ma.pi*Rde*aeff*aeff) * (pp-1.)
nfac = 10.0 
n0 = nfac*n10 ##6.*10**5*Msun/md * 1./(4./3.*ma.pi*(Rrout**3 - Rde**3))


Lav = L0
betst = 0.07#0.06776   ##changed in final version ## gives 0.14 mag amplitdue at alph=1.1
Inc = 0.0#ma.pi/2#ma.acos(0.07/betst)#0.*np.pi/4.
Ombn = 2.*ma.pi/(Rde/c)  ##(2pi/P)
alph = 1.0 # this isfor what might be observed in optical
alphDop = -1.0  #For IR emission this is always -1
Dst = 1.4*10**9*pc2cm




# ### WISE BAND + Observational STUFF
# W1mx = numicron/2.8
# W1mn = numicron/4.0
# W2mx = numicron/3.9
# W2mn = numicron/5.3

# nuVbnd = c/(545.*10**(-7))
# FVbndRel = 3.636*10**(-20)*nuVbnd 
# FW1Rel = 3.09540*10**(-20)*(W1mn + W1mx)/2
# FW2Rel = 1.7187*10**(-20)*(W2mn + W2mx)/2

## Wise band numbers
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3



nuVbnd = c/(5.45*10**(-5))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-21)*8.8560*10**(13)#(W1mn + W1mx)/2
FW2Rel = 1.71787*10**(-21)*6.4451*10**(13)#(W2mn + W2mx)/2



##TABULATE T's and RHSs
print "Creating look up tables"
NT = 15000
RHS_table = np.zeros(NT)
T_table = np.linspace(1., 3000., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)


### PLOT POINTS
Nt=40
tt = np.linspace(0., 2.,       Nt)*2*ma.pi/Ombn


## INTEGRATION LIMTS FOR ALL nu
Nnumn = 0.0
Nnumx = 5.0
numn = Nnumn*numicron
numx = Nnumx*numicron

#numn = W1mn
#numx = W1mx

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
	if (Sphere):
		print "Plottin F_ISO Sphere"
		###########---------------------------###########
		###########-------------Rin-----------###########
		###########---------------------------###########
		if (Plot_R):
			print "Vary R"
			R1 = 1.*Rde
			R2 = 1.*4./5.*Rde
			R3 = 1.*4./3.*Rde
			Omfac = 1.0 
			Ombn = Ombn*Omfac


			tt = np.linspace(0., 2.,       Nt)*2*ma.pi/Ombn

			arg1 = [Lav, Amp, Ombn, t0, n0, R1, pp, thetTst, JJt, aeff, nu0, nne]
			arg2 = [Lav, Amp, Ombn, t0, n0, R2, pp, thetTst, JJt, aeff, nu0, nne]
			arg3 = [Lav, Amp, Ombn, t0, n0, R3, pp, thetTst, JJt, aeff, nu0, nne]

			FI1 = np.zeros(Nt)
			FI2 = np.zeros(Nt)
			FI3 = np.zeros(Nt)

			FsrcI1 = -2.5*np.log10(Fsrc_Iso(tt, Dst, Lav, Amp, Ombn, t0)/FW1Rel)#/FVbndRel)
			#FsrcI1 = Fsrc_Iso(tt, Dst, Lav, Amp, Ombn, t0)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff
								#F_Ring_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
			for i in range(Nt):
				FI1[i] = -2.5*np.log10(F_Sphere_Iso_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
				FI2[i] = -2.5*np.log10(F_Sphere_Iso_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
				FI3[i] = -2.5*np.log10(F_Sphere_Iso_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)
				
				#FI1[i] = F_Sphere_Iso_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff
				#FI2[i] = F_Sphere_Iso_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff
				#FI3[i] = F_Sphere_Iso_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)*(Dst/aeff)**2 * 4.* ma.pi * aeff*aeff


			nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)

			###PLOT###
			plt.figure()
			#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
			#plt.title(r"Sphere, $F^{\rm{src}}_{\rm{iso}}$, $\Omega = %g \times 2\pi c / R_d$" %Omfac)
			#plt.title(r"Sphere, $\Omega = 2\pi c / R_0$")
			plt.title(r"Isotropic, Sphere, $P =  R_0/c$")

			s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=3)

			IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='#1b9e77', linewidth=3)

			IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='#d95f02', linewidth=3)

			IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='#7570b3', linewidth=3)


			plt.grid(b=True, which='both')
			plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{iso}}$', r'$R_d = R_0$',   r'$R_d = 0.8R_0$',   r'$R_d = 1.\bar{33}R_0$'), loc='upper right', fontsize=18)

			plt.xlabel(r"$N_{\rm{orb}}$")
			plt.ylabel("mag")
			plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
			
			plt.ylim(plt.ylim(7.25, 7.65)[::-1])

			plt.tight_layout()

			Savename = "plots/Iso_and_Dop/Sphere/FISO_Sphere_nrm%g_"%nrm+"Om%g_VaryRin_numin%g_numx%g.png" %(Omfac, Nnumn, Nnumx)
			Savename = Savename.replace('.', 'p')
			Savename = Savename.replace('ppng', '.png')
			plt.savefig(Savename)
			#plt.show()
			####-------END Rin-------####



		


	if (Ring):
		print "Plottin F_ISO RINGS"
		###########---------------------------###########
		###########-------------Rin-----------###########
		###########---------------------------###########
		if (Plot_R):
			print "Vary R"
			R1 = Rde
			R2 = 1.25*Rde
			R3 = 1.33*Rde

			arg1 = [Lav, Amp, Ombn, t0, n0, R1, pp, thetTst, JJt, aeff, nu0, nne]
			arg2 = [Lav, Amp, Ombn, t0, n0, R2, pp, thetTst, JJt, aeff, nu0, nne]
			arg3 = [Lav, Amp, Ombn, t0, n0, R3, pp, thetTst, JJt, aeff, nu0, nne]

			FI1 = np.zeros(Nt)
			FI2 = np.zeros(Nt)
			FI3 = np.zeros(Nt)

			FsrcI1 = -2.5*np.log10(Fsrc_Iso(tt, Dst, Lav, Amp, Ombn, t0)/FVbndRel)
								#F_Ring_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
			for i in range(Nt):
				FI1[i] = -2.5*np.log10(F_Ring_Iso_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
				FI2[i] = -2.5*np.log10(F_Ring_Iso_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
				FI3[i] = -2.5*np.log10(F_Ring_Iso_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)

			nrm = np.mean(FsrcI1) - np.mean(FI1)

			###PLOT###
			plt.figure()
			#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
			plt.title(r"Ring, $F^{\rm{src}}_{\rm{iso}}$, $J = %g$ rad" %(JJt))

			s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=3)

			IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='#1b9e77', linewidth=3)

			IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='#d95f02', linewidth=3)

			IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='#7570b3', linewidth=3)


			plt.grid(b=True, which='both')
			plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{iso}}$', r'$R_d = Rd$',   r'$R_d = 1.25Rd$',   r'$R_d = 1.33Rd$'), loc='upper right')

			plt.xlabel(r"$N_{\rm{orb}}$")
			plt.ylabel("mag")
			plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
			#plt.ylim(plt.ylim(11.7, 12.7)[::-1])

			plt.tight_layout()

			Savename = "plots/Iso_and_Dop/Ring/FISO_Ring_nrm%g_"%nrm+"_J%g_VaryRin_numin%g_numx%g.png" %(JJt, Nnumn, Nnumx)
			Savename = Savename.replace('.', 'p')
			Savename = Savename.replace('ppng', '.png')
			plt.savefig(Savename)

			####-------END Rin-------####


		if (Plot_J):
			print "Vary J"
			JJ1 = 0.0
			JJ2 = ma.pi/4
			JJ3 = ma.pi/3
			JJ4 = ma.pi/2

			arg1 = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJ1, aeff, nu0, nne]
			arg2 = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJ2, aeff, nu0, nne]
			arg3 = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJ3, aeff, nu0, nne]
			arg4 = [Lav, Amp, Ombn, t0, n0, Rde, pp, thetTst, JJ4, aeff, nu0, nne]

			FI1 = np.zeros(Nt)
			FI2 = np.zeros(Nt)
			FI3 = np.zeros(Nt)
			FI4 = np.zeros(Nt)

			FsrcI1 = -2.5*np.log10(Fsrc_Iso(tt, Dst, Lav, Amp, Ombn, t0)/FVbndRel)
								#F_Ring_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
			for i in range(Nt):
				FI1[i] = -2.5*np.log10(F_Ring_Iso_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
				FI2[i] = -2.5*np.log10(F_Ring_Iso_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
				FI3[i] = -2.5*np.log10(F_Ring_Iso_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)
				FI4[i] = -2.5*np.log10(F_Ring_Iso_QuadInt(numn, numx, tt[i], Dst, arg4, RHS_table, T_table)/FW1Rel)

			nrm = np.mean(FsrcI1) - np.mean(FI1)

			###PLOT###
			plt.figure()
			#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
			NRd = Rde/3.08e18
			plt.title(r"Ring, $F^{\rm{src}}_{\rm{iso}}$, $R_d = %g$ pc" %(NRd))

			s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=3)

			IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='#1b9e77', linewidth=3)

			IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='#d95f02', linewidth=3)

			IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='#7570b3', linewidth=3)

			IR4 = plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='#e7298a', linewidth=3)


			plt.grid(b=True, which='both')
			plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0], IR4[0]  ], (r'$F^{\rm{src}}_{\rm{iso}}$', r'$J=0$',   r'$J=\pi/4$',   r'$J=\pi/3$',  r'$J=\pi/2$'), loc='upper right')

			plt.xlabel(r"$N_{\rm{orb}}$")
			plt.ylabel("mag")
			plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
			plt.ylim(plt.ylim(11.7, 12.7)[::-1])

			plt.tight_layout()

			Savename = "plots/Iso_and_Dop/Ring/FISO_Ring_nrm%g_"%nrm+"_Rin%g_VaryJ_numin%g_numx%g.png" %(Rde, Nnumn, Nnumx)
			Savename = Savename.replace('.', 'p')
			Savename = Savename.replace('ppng', '.png')
			plt.savefig(Savename)

			####-------END J-------####












if (FDOP):
###########------------------------------------###########
###########				DOPPLER		           ###########
###########------------------------------------###########
	if (Sphere):
		print "Plottin F_DOP Sphere"
		###########---------------------------###########
		###########-------------Rin-----------###########
		###########---------------------------###########
		if (Plot_R):
			print "Vary R"
			R1 = 1.*Rde
			R2 = 1.*4./5.*Rde
			R3 = 1.*4./3.*Rde
			Omfac = 1.0 
			Ombn = Ombn*Omfac

			Inc = 0.0

			tt = np.linspace(0., 2.,       Nt)*2*ma.pi/Ombn


			#Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
			arg1 = [Lav, betst, Inc, Ombn, alphDop, n0, R1, pp, thetTst, JJt, aeff, nu0, nne]
			arg2 = [Lav, betst, Inc, Ombn, alphDop, n0, R2, pp, thetTst, JJt, aeff, nu0, nne]
			arg3 = [Lav, betst, Inc, Ombn, alphDop, n0, R3, pp, thetTst, JJt, aeff, nu0, nne]

			FI1 = np.zeros(Nt)
			FI2 = np.zeros(Nt)
			FI3 = np.zeros(Nt)


			FsrcI1 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)/FW1Rel)#FVbndRel)
								
			for i in range(Nt):
				FI1[i] = -2.5*np.log10(F_Sphere_Dop_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
				FI2[i] = -2.5*np.log10(F_Sphere_Dop_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
				FI3[i] = -2.5*np.log10(F_Sphere_Dop_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)

			nrm = 0.0#np.mean(FsrcI1) - np.mean(FI1)

			###PLOT###
			plt.figure()
			#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
			#plt.title(r"Sphere, $F^{\rm{src}}_{\rm{Dop}}$, $J = %g$ rad" r" $I = %g$ rad" %(JJt, Inc))
			#plt.title(r"Sphere, $\Omega = 2\pi c / R_0$, $I = %g$" %Inc)
			plt.title(r"Doppler, Sphere, $P =  R_0/c$, $I = %g$" %Inc)

			s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=3)

			IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='#1b9e77', linewidth=3)

			IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='#d95f02', linewidth=3)

			IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='#7570b3', linewidth=3)


			plt.grid(b=True, which='both')
			plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{Dop}}$',  r'$R_d = R_0$',   r'$R_d = 0.8R_0$',   r'$R_d = 1.\bar{33}R_0$'), loc='upper right', fontsize=14)

			plt.xlabel(r"$N_{\rm{orb}}$")
			plt.ylabel("mag")
			plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))

			plt.ylim(plt.ylim(7.25, 7.65)[::-1])

			plt.tight_layout()

			Savename = "plots/Iso_and_Dop/Sphere/FDop_R0_Sphere_nrm%g_"%nrm+"_J%g_Inc%g_VaryRin_numin%g_numx%g.png" %(JJt, Inc, Nnumn, Nnumx)
			Savename = Savename.replace('.', 'p')
			Savename = Savename.replace('ppng', '.png')
			plt.savefig(Savename)

			# Savename = "plots/Iso_and_Dop/Sphere/FDop_Sphere_nrm%g_"%nrm+"Om%g_VaryRin_numin%g_numx%g.png" %(Omfac, Nnumn, Nnumx)
			# Savename = Savename.replace('.', 'p')
			# Savename = Savename.replace('ppng', '.png')
			# plt.savefig(Savename)

			####-------END Rin-------####



		###########---------------------------###########
		###########-------------Inc-----------###########
		###########---------------------------###########
		if (Plot_I):
			print "Vary Inc"
			Inc1 = 0.0
			Inc2 = ma.pi/4.
			Inc3 = ma.pi/2.

			#Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
			arg1 = [Lav, betst, Inc1, Ombn, alphDop, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
			arg2 = [Lav, betst, Inc2, Ombn, alphDop, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
			arg3 = [Lav, betst, Inc3, Ombn, alphDop, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]

			FI1 = np.zeros(Nt)
			FI2 = np.zeros(Nt)
			FI3 = np.zeros(Nt)


			#FsrcI1 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FVbndRel)
			#FsrcI2 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc2, Ombn, alph)/FVbndRel)
			#FsrcI3 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc3, Ombn, alph)/FVbndRel)
			

			FsrcI1 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc1, Ombn, alph)/FW1Rel)
			FsrcI2 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc2, Ombn, alph)/FW1Rel)
			FsrcI3 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc3, Ombn, alph)/FW1Rel)

			for i in range(Nt):
									  #F_Ring_Dop_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
				FI1[i] = -2.5*np.log10(F_Sphere_Dop_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
				FI2[i] = -2.5*np.log10(F_Sphere_Dop_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
				FI3[i] = -2.5*np.log10(F_Sphere_Dop_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)

			nrm = np.mean(FsrcI1) - np.mean(FI1)

			###PLOT###
			plt.figure()
			#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
			plt.title(r"Doppler, Sphere, $P = R_0/c$")

			s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='#1b9e77', linewidth=3)
			s2 = plt.plot(tt/(2*np.pi/Ombn), FsrcI2, linestyle = '--', color='#d95f02', linewidth=3)
			s3 = plt.plot(tt/(2*np.pi/Ombn), FsrcI3, linestyle = '--', color='#7570b3', linewidth=3)

			IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='#1b9e77', linewidth=3)

			IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='#d95f02', linewidth=3)

			IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='#7570b3', linewidth=3)


			plt.grid(b=True, which='both')
			plt.legend( [ s1[0],s2[0],s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{Dop}}$', '','', r'$I = 0$',   r'$I = \pi/4$',   r'$I = \pi/2$'), loc='upper right', fontsize=14)

			plt.xlabel(r"$N_{\rm{orb}}$")
			plt.ylabel("mag")
			plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
			plt.ylim(plt.ylim(7.25, 7.65)[::-1])

			plt.tight_layout()

			Savename = "plots/Iso_and_Dop/Sphere/FDop_Sphere_nrm%g_"%nrm+"Rde%g_VaryInc_numin%g_numx%g.png" %(Rde, Nnumn, Nnumx)
			Savename = Savename.replace('.', 'p')
			Savename = Savename.replace('ppng', '.png')
			plt.savefig(Savename)

			####-------END Rin-------####


		if (Plot_alp):
			print "Vary alpha"
			a1 = 2.
			a2 = 3.
			a3 = 4.
			#Omfac = 0.5 
			#Ombn = Ombn*Omfac
			Rde = 1.0*Rde

			Inc = 0.0

			tt = np.linspace(0., 2.,       Nt)*2*ma.pi/Ombn


			#Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
			arg1 = [Lav, betst, Inc, Ombn, a1, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
			arg2 = [Lav, betst, Inc, Ombn, a2, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]
			arg3 = [Lav, betst, Inc, Ombn, a3, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]

			FI1 = np.zeros(Nt)
			FI2 = np.zeros(Nt)
			FI3 = np.zeros(Nt)


			FsrcI1 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, a1)/FW1Rel)#FVbndRel)
			FsrcI2 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, a2)/FW1Rel)#FVbndRel)
			FsrcI3 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, a3)/FW1Rel)#FVbndRel)

								
			for i in range(Nt):
				FI1[i] = -2.5*np.log10(F_Sphere_Dop_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
				FI2[i] = -2.5*np.log10(F_Sphere_Dop_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
				FI3[i] = -2.5*np.log10(F_Sphere_Dop_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)

			nrm1 = np.mean(FsrcI1) - np.mean(FI1)
			nrm2 = np.mean(FsrcI1) - np.mean(FI2)
			nrm3 = np.mean(FsrcI1) - np.mean(FI3)

			###PLOT###
			plt.figure()
			#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
			#plt.title(r"Sphere, $F^{\rm{src}}_{\rm{Dop}}$, $J = %g$ rad" r" $I = %g$ rad" %(JJt, Inc))
			#plt.title(r"Sphere, $\Omega = 2\pi c / R_0$, $I = %g$" %Inc)
			plt.title(r"Doppler, Sphere, $P =  R_0/c$, $I = %g$" %Inc)

			s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='#1b9e77', linewidth=3)
			s2 = plt.plot(tt/(2*np.pi/Ombn), FsrcI2, linestyle = '--', color='#d95f02', linewidth=3)
			s3 = plt.plot(tt/(2*np.pi/Ombn), FsrcI3, linestyle = '--', color='#7570b3', linewidth=3)


			IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm1, color='#1b9e77', linewidth=3)

			IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm2, color='#d95f02', linewidth=3)

			IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm3, color='#7570b3', linewidth=3)

			plt.axhline(y=np.mean(FsrcI1), color='black', linestyle="-")

			plt.grid(b=True, which='both')
			plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{Dop}}$',  r'$\alpha = %g$' %a1,   r'$\alpha = %g$' %a2,   r'$\alpha = %g$' %a3), loc='upper right', fontsize=14)

			plt.xlabel(r"$N_{\rm{orb}}$")
			plt.ylabel("mag")
			plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
			#plt.ylim(plt.ylim(11.7, 12.5)[::-1])

			plt.tight_layout()

			Savename = "plots/Iso_and_Dop/Sphere/FDop_R0_Sphere_nrm%g_"%nrm1+"_J%g_Inc%g_Varyalpha_numin%g_numx%g.png" %(JJt, Inc, Nnumn, Nnumx)
			Savename = Savename.replace('.', 'p')
			Savename = Savename.replace('ppng', '.png')
			plt.savefig(Savename)




	if (Ring):
		print "Plottin F_DOP RINGS"
		###########---------------------------###########
		###########-------------Rin-----------###########
		###########---------------------------###########
		if (Plot_R):
			print "Vary R"
			R1 = Rde
			R2 = 1.25*Rde
			R3 = 1.33*Rde

			#Lavg, bets, incl, Ombin, alphnu, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
			arg1 = [Lav, betst, Inc, Ombn, alphDop, n0, R1, pp, thetTst, JJt, aeff, nu0, nne]
			arg2 = [Lav, betst, Inc, Ombn, alphDop, n0, R2, pp, thetTst, JJt, aeff, nu0, nne]
			arg3 = [Lav, betst, Inc, Ombn, alphDop, n0, R3, pp, thetTst, JJt, aeff, nu0, nne]

			FI1 = np.zeros(Nt)
			FI2 = np.zeros(Nt)
			FI3 = np.zeros(Nt)


			FsrcI1 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)/FVbndRel)
								
			for i in range(Nt):
									  #F_Ring_Dop_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
				FI1[i] = -2.5*np.log10(F_Ring_Dop_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
				FI2[i] = -2.5*np.log10(F_Ring_Dop_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
				FI3[i] = -2.5*np.log10(F_Ring_Dop_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)

			nrm = np.mean(FsrcI1) - np.mean(FI1)

			###PLOT###
			plt.figure()
			#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
			plt.title(r"Ring, $F^{\rm{src}}_{\rm{Dop}}$, $J = %g$ rad" r" $I = %g$ rad" %(JJt, Inc))

			s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=3)

			IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='#1b9e77', linewidth=3)

			IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='#d95f02', linewidth=3)

			IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='#7570b3', linewidth=3)


			plt.grid(b=True, which='both')
			plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0]  ], (r'$F^{\rm{src}}_{\rm{Dop}}$', r'$R_d = Rd$',   r'$R_d = 1.25Rd$',   r'$R_d = 1.33Rd$'), loc='upper right')

			plt.xlabel(r"$N_{\rm{orb}}$")
			plt.ylabel("mag")
			plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
			plt.ylim(plt.ylim(11.7, 12.7)[::-1])

			plt.tight_layout()

			Savename = "plots/Iso_and_Dop/Ring/FDop_Ring_nrm%g_"%nrm+"_J%g_Inc%g_VaryRin_numin%g_numx%g.png" %(JJt, Inc, Nnumn, Nnumx)
			Savename = Savename.replace('.', 'p')
			Savename = Savename.replace('ppng', '.png')
			plt.savefig(Savename)

			####-------END Rin-------####

		if (Plot_J):
			print "Vary J"
			JJ1 = 0.0
			JJ2 = ma.pi/4.
			JJ3 = ma.pi/3.
			JJ4 = ma.pi/2.

			
			arg1 = [Lav, betst, Inc, Ombn, alphDop, n0, Rde, pp, thetTst, JJ1, aeff, nu0, nne]
			arg2 = [Lav, betst, Inc, Ombn, alphDop, n0, Rde, pp, thetTst, JJ2, aeff, nu0, nne]
			arg3 = [Lav, betst, Inc, Ombn, alphDop, n0, Rde, pp, thetTst, JJ3, aeff, nu0, nne]
			arg4 = [Lav, betst, Inc, Ombn, alphDop, n0, Rde, pp, thetTst, JJ4, aeff, nu0, nne]

			FI1 = np.zeros(Nt)
			FI2 = np.zeros(Nt)
			FI3 = np.zeros(Nt)
			FI4 = np.zeros(Nt)

			FsrcI1 = -2.5*np.log10(Fsrc_Dop(tt, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, alph)/FVbndRel)
								#F_Ring_Iso_QuadInt(numin, numax, t, Dist, Aargs, RHStable, Ttable)
			for i in range(Nt):
				FI1[i] = -2.5*np.log10(F_Ring_Dop_QuadInt(numn, numx, tt[i], Dst, arg1, RHS_table, T_table)/FW1Rel)
				FI2[i] = -2.5*np.log10(F_Ring_Dop_QuadInt(numn, numx, tt[i], Dst, arg2, RHS_table, T_table)/FW1Rel)
				FI3[i] = -2.5*np.log10(F_Ring_Dop_QuadInt(numn, numx, tt[i], Dst, arg3, RHS_table, T_table)/FW1Rel)
				FI4[i] = -2.5*np.log10(F_Ring_Dop_QuadInt(numn, numx, tt[i], Dst, arg4, RHS_table, T_table)/FW1Rel)

			nrm = np.mean(FsrcI1) - np.mean(FI1)

			###PLOT###
			plt.figure()
			#plt.title(r"$n_0 = %g n_T$  $J = %g$ rad" %(nfac, JJt))
			NRd = Rde/3.08e18
			plt.title(r"Ring, $F^{\rm{src}}_{\rm{Dop}}$, $R_d = %g$ pc " r" $I = %g$ rad" %(NRd, Inc))

			s1 = plt.plot(tt/(2*np.pi/Ombn), FsrcI1, linestyle = '--', color='blue', linewidth=3)

			IR1 = plt.plot(tt/(2*np.pi/Ombn), FI1+nrm, color='#1b9e77', linewidth=3)

			IR2 = plt.plot(tt/(2*np.pi/Ombn), FI2+nrm, color='#d95f02', linewidth=3)

			IR3 = plt.plot(tt/(2*np.pi/Ombn), FI3+nrm, color='#7570b3', linewidth=3)

			IR4 = plt.plot(tt/(2*np.pi/Ombn), FI4+nrm, color='#e7298a', linewidth=3)


			plt.grid(b=True, which='both')
			plt.legend( [ s1[0], IR1[0], IR2[0], IR3[0], IR4[0]  ], (r'$F^{\rm{src}}_{\rm{Dop}}$', r'$J=0$',   r'$J=\pi/4$',   r'$J=\pi/3$',  r'$J=\pi/2$'), loc='upper right')

			plt.xlabel(r"$N_{\rm{orb}}$")
			plt.ylabel("mag")
			plt.xlim(tt[0]* Ombn/(2.*ma.pi), tt[len(tt)-1] * Ombn/(2.*ma.pi))
			plt.ylim(plt.ylim(11.7, 12.7)[::-1])

			plt.tight_layout()

			Savename = "plots/Iso_and_Dop/Ring/FDop_Ring_nrm%g_"%nrm+"_Rin%g_Inc%g_VaryJ_numin%g_numx%g.png" %(Rde, Inc, Nnumn, Nnumx)
			Savename = Savename.replace('.', 'p')
			Savename = Savename.replace('ppng', '.png')
			plt.savefig(Savename)

			####-------END J-------####

