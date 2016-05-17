import numpy as np
from numpy import *

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt

#import IR_LightEchoes_NewMeth as IRLE
from FluxFuncs_IRLE import *
from ErrFuncs_IRLE import *





###############################
### Define Constants
################################
################################
nne = 1.0
nu0 = numicron/1.5

#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
zPG1302 = 0.2784
L0 = 6.78*10**46 * 1.35
MPGmx = 10**9.4*Msun
RdPG = ma.sqrt(0.1)*2.8 *pc2cm
#OmPG = 2.*ma.pi/(1884.*24.*3600.) * (1.+zPG1302)

OmPG = 2.*ma.pi/(1.87091995e+03*24.*3600.) * (1.+zPG1302)



## TEST VALUES
Lav = L0
betst = 0.1
Inc = ma.acos(0.067/betst)#0.*np.pi/4.
Ombn = OmPG
alph = -3.0

Rde = RdPG
pp = 2.0
thetTst = 0.*np.pi/4
JJt =1.*np.pi/4.
aeff = (c/nu0)/(2.*ma.pi)
#0.16*10**(-4) #(1 micrometer wavelength /(2pi) )


Dst = 1.4*10**9*pc2cm
Rrout = 1.0*Rde



## Wise band numbers
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3

nuVbnd = c/(545*10**(-7))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-20)*(W1mn + W1mx)/2
FW2Rel = 1.7187*10**(-20)*(W2mn + W2mx)/2

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
NT = 10000
RHS_table = np.zeros(NT)
T_table = np.linspace(1., 2000., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)





################################
################################
### PICK TST PARAMS TO PLOT
################################
################################
thetTst = 1.*np.pi/4
JJt =1.*np.pi/2.
sinJJ = ma.sin(JJt)
cosTT = ma.cos(thetTst)
Rin = 1.0 # in units of RdPG

Shell_File = "Testing"
param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$Rin$']
## starting point
#p_tst = [sinJJ,  cosTT,  2.0]
pW1 = [0.05311,     0.90859092,  2.28146583 ]
pW2 = [0.05311,     0.90859092, 10.75429423]


W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 

##opt thin geo thick
#p_thin = [0.7311,  0.80859092,  7.28146583, 1000.]
#p_thin = [6.83121563e-01,   8.88186589e-01,   6.97427899e+00,   1.09843750e+03]
p_thin = [6.47137734e-01,   9.47883340e-01,   6.74388886e+00,   1.23476562e+03]


################################
###############################
### PLOT
################################
################################
Nt=40

ttopt = np.linspace(tsrt[0]-100, t_MJD[len(t_MJD)-1]+100,       Nt)

#Fsrc_Args = [FVbndRel, Lav, betst, OmPG, Dst]
Fsrc_Args = [FVbndRel, Lav, Dst]
Fsrc_p_opt = [ 5.98144879e-02,   6.12468791e-02,   6.55067929e-01,  -3.28334799e-04, 1.87091995e+03]#]1.89037933e+03]

#Fsrc_p_opt  = sc.optimize.fmin(Fsrc_Err2,    fminFsrc_p0, args=(tsrt, Lumsrt, sigLsrt, Fsrc_Args), full_output=1, disp=False)[0]
OmFit = 2.*ma.pi/(Fsrc_p_opt[4]*24.*3600.)* (1.+0.2784)

#opti = -2.5*np.log10(Fsrc((ttopt*3600.*24 - fminFsrc_opt[2]*2.*ma.pi/OmPG), Dst, ma.pi/2., 0.0, fminFsrc_opt[0]*Lav, fminFsrc_opt[1], np.arccos(0.067/fminFsrc_opt[1]), OmPG, 1.1)/FVbndRel)
opti = -2.5*np.log10(Fsrc_Dop( ((ttopt*3600.*24)-Fsrc_p_opt[2]*2.*ma.pi/OmFit), Dst, ma.pi/2., 0.0, Fsrc_p_opt[0]*Lav, Fsrc_p_opt[1], Fsrc_p_opt[3], OmFit, 1.1)/FVbndRel)


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


#W1shell = plt.plot(ttopt, magPoint_OpThick_TorShell(pW1, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
#W2shell = plt.plot(ttopt, magPoint_OpThick_TorShell(pW2, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)


W1shell = plt.plot(ttopt, magPoint_OpThin_TorThick(p_thin, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
W2shell = plt.plot(ttopt, magPoint_OpThin_TorThick(p_thin, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)



plt.grid(b=True, which='both')
		#plt.legend( [ s1[0], IR1[0], s2[0], IR2[0], s3[0], IR3[0]  ], (r'$i=0$','',   r'$i=\pi/4$','',   r'$i=\pi/2$', ''), loc='upper right')

plt.xlabel(r"$t$ [MJD]")
plt.ylabel("mag")
#plt.xlim(52000, 57500)
plt.xlim(3000, max(ttopt))
#plt.ylim(10.5, 11.5)
#plt.ylim(plt.ylim(10.5, 12.3)[::-1])

plt.show()
#plt.savefig("../emcee_data/"+Shell_File+"BestFit.png")
plt.clf()

print "ALL DONE"







