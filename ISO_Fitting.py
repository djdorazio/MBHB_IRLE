import numpy as np
from numpy import *
import scipy as sc
from scipy import *

import emcee

import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
#matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt

from scipy import optimize
from scipy.optimize import fmin

from FluxFuncs_IRLE import *
from ErrFuncs_IRLE import *




################################
###############################
### OPTIONS
################################
################################
Fit_src = True

### Fit a geometrrically (2aeff) thin, optically IR thick model where all emission comes from Rin, don't fit for dust
same_rem = False

## Same as "same_rem" model but fit for dust params = nn, nu0 -> aeff
fit_dust = False
## do the same as fit_dust but allow optically thin to IR
shell_thin = False

## Same as "same_rem" model but allow two diff emission radii for W1 and W2
diff_rem = False

### Fit a geometrically Thick (where tau->1), optically IR thin model with parameters [sinJ, cosT, Rin, n0], don't fit for dust
Opt_Thin = False

same_rem_thin = False



################################
###############################
### Define Constants
################################
################################
nne = 1.
nu0 = numicron/1.5

#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
zPG1302 = 0.2784
L0 = 6.78*10**46 * 1.35
MPGmx = 10**9.4*Msun
RdPG = ma.sqrt(0.1)*2.8 *pc2cm
#OmPG = 2.*ma.pi/(1884.*24.*3600.) * (1.+zPG1302)
#Best fit for source period
OmPG = 2.*ma.pi/(1.87091995e+03*24.*3600.) * (1.+zPG1302)


## TEST VALUES
Lav = L0
betst = 0.068
Inc = ma.acos(0.067/betst)#0.*np.pi/4.
Ombn = OmPG
alph = -2.0

Rde = RdPG
pp = 2.0
thetTst = 1.*np.pi/4
JJt =1.*np.pi/4.
aeff = (c/nu0)/(2.*ma.pi)


Dst = 1.4*10**9*pc2cm
Rrout = 1.0*Rde



## Wise band numbers
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3



nuVbnd = c/(5.45*10**(-5))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-21)*8.8560*10**(13)#(W1mn + W1mx)/2
FW2Rel = 1.71787*10**(-21)*6.4451*10**(13)#(W2mn + W2mx)/2

## PARAMS TO FIT - note these are all params which the optical data does not fit for
#beta0, JJ0, Rin0, nDust0
beta0 = betst
cosJJ0 = ma.cos(JJt)
Rin0 = RdPG









################################
###############################
### IMPORT DATA
################################
################################
print "Importing Data to fit..."

#Import Data to fit
## OPTICAL DATA
Tt   =  (loadtxt("../dat/Lums_PG1302.dat", usecols=[0]))/(1.+zPG1302) 
Lum  =  loadtxt("../dat/Lums_PG1302.dat", usecols=[1])
sigL =  loadtxt("../dat/Lums_PG1302.dat", usecols=[2])

Lum_Mean = mean(Lum)
Lum      = (Lum)# - Lum_Mean)  ## for def of mag - mag0

## put in time order
tLumS   = zip(Tt,Lum,sigL)
tLumS.sort()
TtLumS  = transpose(tLumS)
tsrt    =  TtLumS[0]
Lumsrt  =  TtLumS[1]
sigLsrt =  TtLumS[2]

### IR DATA
t_MJD = (np.genfromtxt("../dat/all.pg1302.txt",usecols=23, comments="|"))/(1.+zPG1302) ## put in binary frame

W1_mag = np.genfromtxt("../dat/all.pg1302.txt",usecols=7, comments="|")
W1_sig = np.genfromtxt("../dat/all.pg1302.txt",usecols=8, comments="|")

W2_mag = np.genfromtxt("../dat/all.pg1302.txt",usecols=10, comments="|")
W2_sig = np.genfromtxt("../dat/all.pg1302.txt",usecols=11, comments="|")


W2_sig = np.genfromtxt("../dat/all.pg1302.txt",usecols=11, comments="|")

##error flags
qual_frame  = np.genfromtxt("../dat/all.pg1302.txt",usecols=19, comments="|")
saa_sep     = np.genfromtxt("../dat/all.pg1302.txt",usecols=21, comments="|")
moon_masked = np.genfromtxt("../dat/all.pg1302.txt",usecols=22, comments="|")

idelq = np.where(qual_frame == 0)[0]
idels = np.where(saa_sep < 0.)[0]
idelm = np.where(moon_masked == 1)[0]

##remove flagged values
W2_mag=np.delete(W2_mag,idelq)
W2_sig=np.delete(W2_sig,idelq)
t_MJD = np.delete(t_MJD,idelq)
W1_mag=np.delete(W1_mag,idelq)
W1_sig=np.delete(W1_sig,idelq)

W2_mag=np.delete(W2_mag,idels)
W2_sig=np.delete(W2_sig,idels)
t_MJD = np.delete(t_MJD,idels)
W1_mag=np.delete(W1_mag,idels)
W1_sig=np.delete(W1_sig,idels)

W2_mag=np.delete(W2_mag,idelm)
W2_sig=np.delete(W2_sig,idelm)
t_MJD = np.delete(t_MJD,idelm)
W1_mag=np.delete(W1_mag,idelm)
W1_sig=np.delete(W1_sig,idelm)



###### get average value for each cluser of data points in time
iseg = []
iseg.append(-1)
for i in range(len(t_MJD)-1):
	if (abs((t_MJD[i] - t_MJD[i+1])) > 10 ):
		iseg.append(i)
iseg.append(len(t_MJD)-1)

t_avg = []
W1_avg = []
W2_avg = []
W1_avsg = []
W2_avsg = []


for i in range(0 , len(iseg)-1):
	t_avg.append(np.mean(t_MJD[iseg[i]+1:iseg[i+1]]))

	W1_avg.append(np.mean(W1_mag[iseg[i]+1:iseg[i+1]]))
	W2_avg.append(np.mean(W2_mag[iseg[i]+1:iseg[i+1]]))

	#W1_avsg.append((max(W1_mag[iseg[i]+1:iseg[i+1]]) - min(W1_mag[iseg[i]+1:iseg[i+1]]))/6.)
	#W2_avsg.append((max(W2_mag[iseg[i]+1:iseg[i+1]]) - min(W2_mag[iseg[i]+1:iseg[i+1]]))/6.)
	Nseg = len(W1_sig[iseg[i]+1:iseg[i+1]])
	W1_avsg.append(np.sqrt(sum( (W1_sig[iseg[i]+1:iseg[i+1]])**2 ))/Nseg)
	W2_avsg.append(np.sqrt(sum( (W2_sig[iseg[i]+1:iseg[i+1]])**2 ))/Nseg)

## APPEND THE ONE AKARI DATA POINT
#MJD 51537.0 (ISO), 55022.5 (Akari)
#t_avg.append(51537.0-50000.0)
t_avg.append(55022.5/(1.+zPG1302))
#W1 NaN, 11.3008
W1_avg.append(11.3008)
#W1_error NaN, 0.0603544
W1_avsg.append(0.0603544)
#W2 10.3370, 10.2411
W2_avg.append(10.2411)
#W2_error 0.217200, 0.0691765
W2_avsg.append(0.0691765)

t_avg = np.array(t_avg)
W1_avg = np.array(W1_avg)
W2_avg = np.array(W2_avg)
W1_avsg = np.array(W1_avsg)
W2_avsg = np.array(W2_avsg)





################################
###############################
### Temp LOOKUP TABLE
################################
################################
#if (ShellFit or ThickFit or fit_both or pltShell or pltThick or pltboth):
#Set up look up tables
##TABULATE T's and RHSs
print "Creating Temp look up tables..."
NT = 2000
RHS_table = np.zeros(NT)
T_table = np.linspace(1., 2000., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)



################################
###############################
### FMIN FITTING
################################
################################
###Fit optical
if (Fit_src):
	Shell_File = "ISO_FsrcOptical_Lfac_Amp_phs_Pday_"
	param_names = ['Lfac','Amp', 'phase', 'Prd']
	## starting point
	Fsrc_ISO_p0 = [0.0597279747, 0.139181205, 0.688098413, 1871.99573]
	## args of non chanigng parameters to pass
	Fsrc_Args = [FVbndRel, Lav, Dst]
	## optimize with fmin
	Fsrc_p_opt  = sc.optimize.fmin(Fsrc_ISO_Err2,    Fsrc_ISO_p0, args=(tsrt, Lumsrt, sigLsrt, Fsrc_Args), full_output=1, disp=False)[0]

	pW1 = Fsrc_p_opt
	pW2 = Fsrc_p_opt
	ps = Fsrc_p_opt



if (fit_dust):
	sinJJ = ma.sin(JJt)
	cosTT = ma.cos(thetTst)
	Rin = 1.0 # in units of RdPG

	if (shell_thin):
		Shell_File = "ISO_OptIRTHIN_Dust_Fit_TorusShell_J_ThT_Rin_nn_n0_"
		param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$Rin$',r'Amp', r'$k$' , r'$nu_0$']
		## starting point
		#OpThick_TorShell_p0 = [sinJJ, cosTT, Rin]
		#OpThin_TorShell_dust_p0 = [-5.81648541e-01,   0.7,   9.28085681e+00,   3.62563978e-01, 4.12632801e+14]
		#chi2 278 below
		OpThin_TorShell_dust_p0 =  [-6.60098675e-01,   4.57835525e-01,   8.67955344e+00,   0.36, 5.12336423e-01, 3.35224080e+14]
		## args of non chanigng parameters to pass
		Fsrc_ISO_p0 = [0.0597279747, 0.139181205, 0.688098413, 1871.99573]
		Ombn =	0.0#2.*ma.pi/(Fsrc_ISO_p0[3]*24.*3600.) * (1.+0.2784)
		t0   = 0.0#Fsrc_ISO_p0[2] * 2.*ma.pi/Ombn
		#Amp  = Fsrc_ISO_p0[1]

		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, pp, aeff,    t0] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, pp, aeff,    t0] 
		## optimize with fmin
		OpThin_TorShell_p_opt  = sc.optimize.fmin(ISO_OpThin_TorShell_Err2,    OpThin_TorShell_dust_p0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.1)[0]
		pW1 = OpThin_TorShell_p_opt
		pW2 = OpThin_TorShell_p_opt
		ps = OpThin_TorShell_p_opt


	else:
		Shell_File = "ISO_OptIRTHICK_Dust_Fit_TorusShell_J_ThT_Rin_nn_n0_"
		param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$Rin$', r'$nn$' , r'$nu_0$']
		## starting point
		#OpThick_TorShell_p0 = [sinJJ, cosTT, Rin]
		#OpThick_TorShell_dust_p0 = [-5.81648541e-01,   0.7,   9.28085681e+00,   3.62563978e-01, 4.12632801e+14]
		OpThick_TorShell_dust_p0 = [0.7,   6.25236770e-01,   6.60443550e+00,   0.0, 4.50861670e+14]
		## args of non chanigng parameters to pass
		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  betst] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  betst] 
		## optimize with fmin
		OpThick_TorShell_p_opt  = sc.optimize.fmin(OpThick_TorShell_dustP_Err2,    OpThick_TorShell_dust_p0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.1)[0]
		pW1 = OpThick_TorShell_p_opt
		pW2 = OpThick_TorShell_p_opt
		ps = OpThick_TorShell_p_opt

if (same_rem_thin):
	sinJJ = ma.sin(JJt)
	cosTT = ma.cos(thetTst)
	Rin = 1.0 # in units of RdPG

	Shell_File = "ISO_OptThin_TorusShell_J_ThT_Rin_Amp_"
	param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$Rin$', 'Amp']
	#OpThin_TorShell_p0 =  [-6.60098675e-01,   4.57835525e-01,   8.67955344e+00,   0.36]

	
	#chi2 = 242
	OpThin_TorShell_p0 =  [-0.94789069,   0.91328937,  10.23015572,   0.09667759]


	## args of non chanigng parameters to pass
	Fsrc_ISO_p0 = [0.0597279747, 0.139181205, 0.688098413, 1871.99573]
	Ombn =	0.0#2.*ma.pi/(Fsrc_ISO_p0[3]*24.*3600.) * (1.+0.2784)
	t0   = 0.0#Fsrc_ISO_p0[2] * 2.*ma.pi/Ombn
	#Amp  = Fsrc_ISO_p0[1]

	W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, pp, aeff, nu0, nne,   t0] 
	W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, pp, aeff, nu0, nne,   t0] 
	## optimize with fmin
	OpThin_TorShell_p_opt  = sc.optimize.fmin(ISO_OpThin_TorShell_Err2,    OpThin_TorShell_p0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.1)[0]
	pW1 = OpThin_TorShell_p_opt
	pW2 = OpThin_TorShell_p_opt
	ps = OpThin_TorShell_p_opt




######################################
### FIT FOR W1 and W2 at different Rin
######################################
if (diff_rem):
	Shell_File = "ISO_OptThin_TorusShell_J_ThT_Rin1_Rin2_"
	param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$Rin1$',r'$Rin2$']
	DiffR_OpThick_TorShell_p0 = [-0.06475812,   0.90372837,   2.39385413,  10.32260737]
	#ch2 = 249.227189075
	## args of non chanigng parameters to pass
	W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
	W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
	## optimize with fmin
	DiffR_OpThick_TorShell_p_opt  = sc.optimize.fmin(DiffR_OpThick_TorShell_Err2,    DiffR_OpThick_TorShell_p0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.1)[0]
	pW1 = [DiffR_OpThick_TorShell_p_opt[0], DiffR_OpThick_TorShell_p_opt[1], DiffR_OpThick_TorShell_p_opt[2]]
	pW2 = [DiffR_OpThick_TorShell_p_opt[0], DiffR_OpThick_TorShell_p_opt[1], DiffR_OpThick_TorShell_p_opt[3]]
	ps = DiffR_OpThick_TorShell_p_opt



if (Opt_Thin):
	Shell_File = "ISO_OptThn_TorusThick_J_ThT_Rin_n0_"
	param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$R_{in}$', r'$n_0$', r'$A$']
	## starting point
	#182.3 below
	#OpThin_TorThick_p0 = [6.72054830e-01,   9.59726315e-01,   6.84491017e+00,   1.24892755e+03]
	# yeti, chi2 = 
	#ISO_OpThin_TorThick_p0 = [6.79489126e-01,   9.84816675e-01,   6.78328479e+00,   1.20359057e+03, 0.35]
	#chi2 174 -[6.71240296e-01   9.68638091e-01   6.84514215e+00   1.26185273e+03, 3.83701459e-01]
	#chi2 = 140
	#ISO_OpThin_TorThick_p0 = [6.63012874e-01,   9.52632866e-01,   6.93571877e+00,   1.25865725e+03, 3.47533200e-01]		
	#chi2 140
	#ISO_OpThin_TorThick_p0 = [6.84661480e-01,   9.43611434e-01,   6.96853459e+00,   1.21721076e+03, 3.44821368e-01]
	#chi2 122
	ISO_OpThin_TorThick_p0 = [6.81276630e-01,   9.39613783e-01,   6.94579449e+00,   1.23506837e+03, 3.43593216e-01]


	Fsrc_ISO_p0 = [0.0597279747, 0.139181205, 0.688098413, 1871.99573]
	Ombn =	0.0#2.*ma.pi/(Fsrc_ISO_p0[3]*24.*3600.) * (1.+0.2784)
	t0   = 0.0#Fsrc_ISO_p0[2] * 2.*ma.pi/Ombn
	#Amp  = Fsrc_ISO_p0[1]

	W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, pp, aeff, nu0, nne,    t0] 
	W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, pp, aeff, nu0, nne,    t0] 
	## optimize with fmin
	ISO_OpThin_TorThick_p_opt  = sc.optimize.fmin(ISO_OpThin_TorThick_Err2,    ISO_OpThin_TorThick_p0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.1)[0]
	pW1 = ISO_OpThin_TorThick_p_opt
	pW2 = ISO_OpThin_TorThick_p_opt
	ps = ISO_OpThin_TorThick_p_opt





###############################
### Write FMIN VALUES OT FILE
################################
filename = "fmin"+Shell_File+"_results.txt"
print "Printing Results"
target = open(filename, 'w')
#target.truncate()

	
for i,name in enumerate(param_names):
	target.write("W1: {name}: {0:.4f}".format(ps[i], name=name))
	target.write("\n")

target.write("\n")

for i,name in enumerate(param_names):
	target.write("W2: {name}: {0:.4f}".format(ps[i], name=name))
	target.write("\n")
	
		

target.close

filename = "fmin"+Shell_File+"_results.txt"
print "Printing Results"
target = open(filename, 'w')
#target.truncate()

	
	
for i,name in enumerate(param_names):
	target.write("W1: {name}: {0:.4f}".format(ps[i], name=name))
	target.write("\n")

target.write("\n")

for i,name in enumerate(param_names):
	target.write("W2: {name}: {0:.4f}".format(ps[i], name=name))
	target.write("\n")

	
		

target.close



################################
###############################
### PLOT
################################
################################
Nt=40

ttopt = np.linspace(tsrt[0]-100, t_MJD[len(t_MJD)-1]+100,       Nt)


Fsrc_p_opt = [0.0597279747, 0.139181205, 0.688098413, 1871.99573]
Ombn =	2.*ma.pi/(Fsrc_p_opt[3]*24.*3600.) * (1.+0.2784)
t0   = Fsrc_p_opt[2] * 2.*ma.pi/Ombn
opti = -2.5*np.log10(Fsrc_Iso((ttopt*3600.*24.-t0), Dst, Fsrc_p_opt[0]*Lav, Fsrc_p_opt[1], Ombn, t0)/FVbndRel) 

ttopt = (ttopt*(1.+zPG1302) - 50000)
t_avg = (t_avg*(1.+zPG1302) - 50000)
tsrt  = (tsrt*(1.+zPG1302) - 50000)
t_MJD = (t_MJD*(1.+zPG1302) - 50000)


print "PLOTTING"

plt.figure()
plt.errorbar(tsrt, Lumsrt-3.0, yerr=sigL, linestyle="none", color = "blue", alpha=0.5) #alpha=0.1

Fs = plt.plot(ttopt, opti-3.0, linestyle = '--', color='blue', linewidth=2)

W1dat   = plt.errorbar(t_MJD, W1_mag, yerr=W1_sig, linestyle="none", color='orange', alpha=1., elinewidth=1.5)
W2dat   = plt.errorbar(t_MJD, W2_mag+0.5, yerr=W2_sig, linestyle="none", color='red', alpha=1., elinewidth=1.5)

W1av   = plt.errorbar(t_avg, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)
W2av   = plt.errorbar(t_avg, W2_avg+0.5, yerr=W2_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)


if (diff_rem or same_rem):
	W1shell = plt.plot(ttopt, magPoint_OpThick_TorShell(pW1, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
	W2shell = plt.plot(ttopt, magPoint_OpThick_TorShell(pW2, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
if (Opt_Thin):
	W1shell = plt.plot(ttopt, ISO_magPoint_OpThin_TorThick(pW1, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
	W2shell = plt.plot(ttopt, ISO_magPoint_OpThin_TorThick(pW2, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
if (fit_dust):
	if (shell_thin):
		W1shell = plt.plot(ttopt, ISO_magPoint_OpThin_TorThin(pW1, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, ISO_magPoint_OpThin_TorThin(pW2, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
	else:
		W1shell = plt.plot(ttopt, magPoint_OpThick_TorShell_dustP(pW1, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_OpThick_TorShell_dustP(pW2, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
if (same_rem_thin):
	W1shell = plt.plot(ttopt, ISO_magPoint_OpThin_TorShell(pW1, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
	W2shell = plt.plot(ttopt, ISO_magPoint_OpThin_TorShell(pW2, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)



plt.grid(b=True, which='both')
		#plt.legend( [ s1[0], IR1[0], s2[0], IR2[0], s3[0], IR3[0]  ], (r'$i=0$','',   r'$i=\pi/4$','',   r'$i=\pi/2$', ''), loc='upper right')

plt.xlabel(r"$t$ [MJD]")
plt.ylabel("mag")
#plt.xlim(52000, 57500)
plt.xlim(3000, max(ttopt))
#plt.ylim(10.5, 11.5)
plt.ylim(plt.ylim(10.5, 12.3)[::-1])

		#plt.show()
plt.savefig("../emcee_data/"+Shell_File+"BestFit.png")
plt.clf()

print "ALL DONE"







