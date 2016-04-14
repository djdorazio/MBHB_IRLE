#Import python modules
import cPickle as pickle
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

from IR_LightEchoes_NewMeth import *

###OPTIONS

## JUST PLOT DONT FIT
NoFit = False   ## jsut plot - turn on one of the plt... options below
pltShell = False
rem_is_Rin = False ### fit for rem and Rin in shell model??
#
pltboth  = False
pltThick = False



## WHAT KIND OF FITTING? ALL TURNED OFF IF NOFIT SET ABOVE
emcee_Fit = False # use emcee to fit
## multiprocessing
NThread = 48
#
fmin_Fit = True   # use simple fmin to fit

W1fit = False  ### fit only W1fit, Thick or Shell
W2fit = False  ### fit only W1fit, Thick or Shell
fit_both = True  ### fit using chi2 for both W1 and W2 - SHELL ONLY - 4 or 6 params depends on rem_is_Rin
				### if rem_is_Rin is FALSE, then 6 params [sinJ, cosT, rem1, rem2, Rin, n0] where rems are
				### where emission for W1 and W2 comes from and Rin is the sublimation radius rem>Rin
if (fit_both):
	W1fit = False
	W2fit = False


Fit_Src = False     ## fmin fit for Lfrac, beta, phase, and Inc to fit optical data
sinFit_Src = False  ## emcee fit a sin curve to optical data
SinFit = False      ## emcee fit a sin curve to IR data (W1 and W2)
No_Prd = True       ## fix the period at 1884/(1+z) ~ 1474 for fitting


ShellFit = False #-> make false if fitboth is on - this fits for W1 and W2 independnetly and in addition to fitboth
ThickFit = False # Fits W1 or W2 as selected above for THick model

mpi_it = False  #defunct


if (NoFit):
	Shell_File = "NoFit"
	emcee_Fit = False
	fmin_Fit = True

	W1fit = False
	W2fit = False
	fit_both = False

	Fit_Src = False
	sinFit_Src = False
	ShellFit = False
	ThickFit = False
	SinFit   = False
	No_Prd = True


#Define Constants
nne = 1.
nu0 = numicron

#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
zPG1302 = 0.2784
Omb = 1./(1*yr2sec)
L0 = 6.78*10**46 * 1.35
MPGmx = 10**9.4*Msun
Ryr = c*yr2sec
RdPG = ma.sqrt(0.1)*2.8 *pc2cm
OmPG = 2.*ma.pi/(1884.*24.*3600.) * (1.+zPG1302)
#alphnu = 0.0

Rorb = c*2.*np.pi/Omb
Ompc = 2.*np.pi*c/pc2cm/2.


## TEST VALUES
Lav = L0
betst = 0.08
Inc = ma.acos(0.067/betst)#0.*np.pi/4.
Ombn = OmPG
alph = -2.0

Rde = RdPG
pp = 2.0
thetTst = 1.*np.pi/4
JJt =4.*np.pi/8
aeff = 0.16*10**(-4) #(0.1 micrometer is an average ISM dust grain size)


Dst = 1.4*10**9*pc2cm
Rrout = 100.0*Rde

md = 10**(-14)
Mdust = 6.*10**5


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
#nDust0 = Mdust*Msun/md * 1./(4./3.*ma.pi*(Rrout**3 - Rde**3))
n0 = 10.0/(ma.pi*Rde*aeff*aeff) * (pp-1.)


# p0 = [Amg, Prd, phs, mag0]
Sinp0_W1 = [0.0887, 1884., 0.01, 11.3]
Sinp0_W2 = [0.1, 1884., 0.01, 10.3]
#no Period fit
if (No_Prd):
	SinPrd = 1884./(1+zPG1302)  ##compare in binary frame
	Sinp0_W1 = [0.0887, 347.7479, 11.3914]
	Sinp0_W2 = [0.0795, 399.3936, 10.3067]



if (ShellFit):
	if (rem_is_Rin):
		#p0 = [sinJ, costheta_T, Rin, n0]
		ShW1_p0_0  = [0.0010,   0.7437, 1.4572,  0.5496] 
		ShW2_p0_0  = [0.0008,   0.6142, 1.2751,  2.9455]
		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
	else:
		#p0 = [sinJ, costheta_T, rem, Rin, n0]
		ShW1_p0_0  = [0.0010,   0.7437, 1.0,   1.0,  0.5496] 
		ShW2_p0_0  = [0.0008,   0.6142, 0.01,  1.0, 2.9455]
		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
if (fit_both):
	if (rem_is_Rin):
		#p0 = [sinJ, costheta_T, R1, n0]
		Shboth_p0_0 = [0.0004,   0.6356, 2.5976,  0.4841] # fit both withh all same parameters (from same inner edge)	
		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst]
	else:
		Shboth_p0_0 = [0.0004,   0.6356, 1.0, 0.01, 1.0,  0.4841]
		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst]
		#ShW1_p0_0 = [Shboth_p0_0[0], Shboth_p0_0[1], Shboth_p0_0[2], Shboth_p0_0[4], Shboth_p0_0[5]]
		#ShW2_p0_0 = [Shboth_p0_0[0], Shboth_p0_0[1], Shboth_p0_0[3], Shboth_p0_0[4], Shboth_p0_0[5]]
if (ThickFit):
	#p0 = [cosJ, costheta_T, Rin, p, n0]
	#ShW1_p0_0  = [ 0.0016,  0.7, 2.0,  1.0]
	#ShW2_p0_0  = [ 0.0016,  0.7, 2.0,  1.0]
	ShW1_p0_0  = [ 0.0016,  0.7, 1.0, 2.0,  10.0]
	ShW2_p0_0  = [ 0.0016,  0.7, 1.0, 2.0,  10.0]
	W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, Rrout,  aeff, nu0, nne, betst] 
	W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, Rrout,  aeff, nu0, nne, betst] 
if (NoFit):
	#ShW1_p0_0  = [ 0.001,   0.6905, 1.4392,  0.5880]
	#ShW2_p0_0  = [ 0.0009,  0.6035, 1.0947,  2.5117]
	ShW1_p0_0  = [ 0.0016,  0.7, 2.0,  1.0]
	ShW2_p0_0  = [ 0.0016,  0.7, 2.0,  1.0]
	if (pltboth): 
		Shell_File = "ShellBoth_SameR_xefix"
		p1both = [0.0004,   0.6356, 2.5976,  0.4841] # fit both withh all same parameters (from same inner edge)
		p2both = p1both
		ShW1_p0_0 = p1both
		ShW2_p0_0 = p1both
		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
	if (pltShell):
		if (rem_is_Rin):
			Shell_File = "Shell_rem_is_Rin_xefix"
			ShW1_p0_0  = [0.0004,   0.6356, 2.5976,  0.4841]
			ShW2_p0_0  = [0.0004,   0.6356, 2.5976,  0.4841]
			W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout, aeff, nu0, nne, betst] 
			W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout, aeff, nu0, nne, betst]
		else:
			Shell_File = "Shell_rem_xefix"
			#p0 = [sinJ, CosT, rem1, rem2, Rin, n0]
			SHboth_p0  = [0.0016,   0.6356, 1.0, 0.01, 1.0,  10.0]
			ShW1_p0_0  = [SHboth_p0[0], SHboth_p0[1],SHboth_p0[2],SHboth_p0[4], SHboth_p0[5]]
			ShW2_p0_0  = [SHboth_p0[0], SHboth_p0[1],SHboth_p0[3],SHboth_p0[4], SHboth_p0[5]]
			W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout, aeff, nu0, nne, betst] 
			W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout, aeff, nu0, nne, betst]

	if (pltThick):
		Shell_File = "Thick_xefix"
		#ShW1_p0_0  = [ 1.0,  0.8957, 0.1, 0.5254,  0.1163]
		#ShW2_p0_0  = [ 1.0,  0.8957, 0.1, 0.5254,  0.1163]
		ShW1_p0_0  = [ 0.0011,  0.8164, 1.9627, 0.6721,  10.0]
		ShW2_p0_0  = [ 0.0011,  0.8164, 1.9627, 0.6721,  10.0]
		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, Rrout, aeff, nu0, nne, betst] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, Rrout, aeff, nu0, nne, betst] 

#Targs = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]






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

#import sys
#sys.exit(0)
### averaging data ^####

if (ShellFit or ThickFit or fit_both or pltShell or pltThick or pltboth):
	#Set up look up tables
	##TABULATE T's and RHSs
	print "Creating Temp look up tables..."
	NT = 10000
	RHS_table = np.zeros(NT)
	T_table = np.linspace(1., 2000., NT)
	for i in range(NT):
		RHS_table[i] = T_RHS(T_table[i], nu0, nne)


##-------- more advanced error fitting ----------##
# ### MCMC - set up Errfuncs
# def Cov_All(t, Lumsrt, sigLsrt, sigG, tauD):
# 	sigG2 = sigG*sigG
# 	Cov_DRW = zeros([len(Lumsrt), len(Lumsrt)])
# 	Cov_phot = zeros([len(Lumsrt), len(Lumsrt)])
# 	for i in range (0, len(Lumsrt)):
# 		for j in range(0, len(Lumsrt)):
# 			Cov_DRW[i][j] = sigG2 * exp(-abs(t[i] - t[j])/(tauD*day*(1.+0.2784)))  ## EQN 4 Zu, Kochanek , Peterson 2011 - convert to obs frame (1+z)
# 			if (i==j): 
# 				Cov_phot[i][j] = sigLsrt[i]*sigLsrt[j]
# 			else:	
# 				Cov_phot[i][j] = 0.0	
				
	
# 	CovT = Cov_DRW + Cov_phot
# 	from scipy import linalg
# 	CovT_I = linalg.inv(CovT)
# 	return {'CovT_I':CovT_I, 'CovT':CovT, 'Cov_DRW':Cov_DRW, 'Cov_phot':Cov_phot}


# def CovErr2(params, t, THEargs, RHStable, Ttable, y, dy):
# 	#print "EVAL", p
# 	#print "EVAL"
# 	sigD = p[len(p)-2]
# 	tauD = p[len(p)-1]
# 	#if (sigD>sigDMAX or tauD < tauDmin or tauD > tauDMAX):
# 	#comment below only for emcee!!!!!!!!!
# 	#if (sigD<0. or tauD < 0.):
# 	#	LnP = -np.inf
# 	#else:
# 	YY = (y - magPoint(params, t, THEargs, RHStable, Ttable))
# 	YY.shape = (len(y),1) #Make column vector
# 	YYT=YY.T
		
# 	CovAll = Cov_All(t, y, dy, sigD, tauD)
# 	NN = np.mat(CovAll['Cov_phot'])
# 	SS = np.mat(CovAll['Cov_DRW'])
		
# 	SNeig = linalg.eig(SS*NN)
# 	SNlogDet = np.sum(np.log(SNeig[0]))
# 	chi2 = YYT.dot(CovT_Inv(t, y, dy, sigD, tauD)).dot(YY)
		
# 	LnP = -0.5*SNlogDet - 0.5*chi2
		
# 	#print(chi2)
# 	#print LnP
# 	return -LnP
##-------- more advanced error fittingabove  ----------##
import time
print "start timing"


### best fit fits both W1 and W2
def Shell_RegErr2(p, t, THEargs, RHStable, Ttable, y, dy,rem_is_Rin):
	print "EVAL", p
	t1=time.clock()
	chi = (y - magPoint_Shell(p, t, THEargs, RHStable, Ttable, rem_is_Rin)) / dy
	#nLnP = sum(chi*chi)
	t2=time.clock()
	#print(chi2)
	print(t2-t1)
	return sum(chi*chi)

	
def ShellBoth_RegErr2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2, rem_is_Rin):
	print "EVAL", p
	t1=time.clock()
	#p0 = [cosJ, cosT, rem1, rem2, Rin, n0]
	if (rem_is_Rin):
		chi1 = (y1 - magPoint_Shell(p, t, THEargs1, RHStable, Ttable, rem_is_Rin)) / dy1
		chi2 = (y2 - magPoint_Shell(p, t, THEargs2, RHStable, Ttable, rem_is_Rin)) / dy2
	else:
		p1 = [p[0], p[1], p[2], p[4], p[5]]
		p2 = [p[0], p[1], p[3], p[4], p[5]]
		chi1 = (y1 - magPoint_Shell(p1, t, THEargs1, RHStable, Ttable, rem_is_Rin)) / dy1
		chi2 = (y2 - magPoint_Shell(p2, t, THEargs2, RHStable, Ttable, rem_is_Rin)) / dy2
	#nLnP = sum(chi*chi)
	t2=time.clock()
	#print(chi2)
	print(t2-t1)
	return sum(chi1*chi1) + sum(chi2*chi2)


def Thick_RegErr2(p, t, THEargs, RHStable, Ttable, y, dy):
	print "EVAL", p
	#t1=time.clock()
	chi = (y - magPoint_Thick(p, t, THEargs, RHStable, Ttable)) / dy
	#nLnP = sum(chi*chi)
	#t2=time.clock()
	#print(chi2)
	#print(t2-t1)
	return sum(chi*chi)

def Thick_RegErr2_fmin(p, t, THEargs, RHStable, Ttable, y, dy):
	print "EVAL", p
	#t1=time.clock()
	chi = (y - magPoint_Thick_fmin(p, t, THEargs, RHStable, Ttable)) / dy
	#nLnP = sum(chi*chi)
	#t2=time.clock()
	#print(chi2)
	#print(t2-t1)
	return sum(chi*chi)



def sinPoint(params, t):
	if (No_Prd):
		Prd = SinPrd
		Amp, phs, mag0 = params
	else:
		Amp, Prd, phs, mag0 = params
	#Amp, phs, mag0 = params
	#Prd=1884.
	return Amp*np.sin( 2.*ma.pi/Prd*(t - phs) - np.pi/2) + mag0 


def SinErr2(p, t, y, dy):
	print "EVAL", p
	chi = (y - sinPoint(p, t) )/ dy
	return sum(chi*chi)
	#print(chi2)
	#return nLnP


def Fsrc_Err2(p, t, y, dy):
	print "EVAL", p
	Lfac, bets, phs, incl = p
	#incl = np.arccos(0.067/bets)
	alphnu =1.1
	Ombin =	OmPG
	chi = (y - -2.5*np.log10(Fsrc((t-phs*2.*ma.pi/Ombin), Dst, ma.pi/2., 0.0, Lfac*Lav, bets, incl, Ombin, alphnu)/FVbndRel) )/ dy
	return sum(chi*chi)
	#print(chi2)
	#return nLnP



### MCMC - Set up priors
def ln_prior(params, rem_is_Rin):
			#beta, cosJJ, Rin, thetT, n0 = p
			if (rem_is_Rin):
				sinJJ, cosTT, Rin, n0 = params
					
				if sinJJ < -1 or sinJJ > 1:
					return -np.inf

				if cosTT < 0 or cosTT > 1:
					return -np.inf
						
				if Rin <= 0.0:
					return -np.inf

				if n0 <= 0.0:
					return -np.inf
			else:
				sinJJ, cosTT, rem1, rem2, Rin, n0 = params
					
				if sinJJ < -1 or sinJJ > 1:
					return -np.inf

				if cosTT < 0 or cosTT > 1:
					return -np.inf

				if rem1 <= 0.0:
					return -np.inf

				if rem2 <= 0.0:
					return -np.inf
						
				if Rin <= 0.0:
					return -np.inf

				if n0 <= 0.0:
					return -np.inf
					
				
			return 0.

def ln_THprior(params):
			sinJJ, cosTT, Rin, pp, n0 = params
					
			if sinJJ < -1 or sinJJ > 1:
				return -np.inf

			if cosTT < 0 or cosTT > 1:
				return -np.inf
					
			if Rin <= 0.0:
				return -np.inf

			if pp <= 0.0:
				return -np.inf

			if n0 <= 0.0:
				return -np.inf
					
				
			return 0.

def ln_Sinprior(p):
	if (No_Prd):
		Prd = SinPrd
		Amp, phs, mag0 = p
		if Amp < 0.:
			return -np.inf

		if phs < 0.0 or phs > Prd:
			return -np.inf
	else:
		Amp, Prd, phs, mag0 = p	
		if Amp < 0.:
			return -np.inf
					
		if Prd < 0.:
			return -np.inf

		if phs < 0.0 or phs > Prd:
			return -np.inf
					
				
	return 0.

### MCMC - Set up posteriors
def ln_Shlikelihood(p, t, Wargs, RHStable, Ttable, y, dy, rem_is_Rin):
			return -(Shell_RegErr2(p, t, Wargs, RHStable, Ttable, y, dy, rem_is_Rin))


def ln_ShBothlikelihood(p, t, W1args, W2args, RHStable, Ttable, y1, dy1,y2, dy2, rem_is_Rin):
			return -(ShellBoth_RegErr2(p, t, W1args, W2args, RHStable, Ttable, y1, dy1,y2, dy2, rem_is_Rin))


def ln_Thlikelihood(p, t, Wargs, RHStable, Ttable, y, dy):
			return -(Thick_RegErr2(p, t, Wargs, RHStable, Ttable, y, dy)) #+ RegErr2(p, t, W2args, RHStable, Ttable, y2, dy2))



def ln_Shposterior(p, t, Wargs, RHStable, Ttable, y, dy, rem_is_Rin):
			ln_p = ln_prior(p, rem_is_Rin)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_Shlikelihood(p, t, Wargs, RHStable, Ttable, y, dy, rem_is_Rin)
			return ln_l + ln_p


def ln_ShBothPosterior(p, t, W1args, W2args, RHStable, Ttable, y1, dy1,y2, dy2, rem_is_Rin):
			ln_p = ln_prior(p, rem_is_Rin)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_ShBothlikelihood(p, t, W1args, W2args, RHStable, Ttable, y1, dy1,y2, dy2, rem_is_Rin)
			return ln_l + ln_p


def ln_Thposterior(p, t, Wargs, RHStable, Ttable, y, dy):
			ln_p = ln_Thprior(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_Thlikelihood(p, t, Wargs, RHStable, Ttable, y, dy)
			return ln_l + ln_p


def ln_Sinlikelihood(p, t, y, dy):
			return -(SinErr2(p, t, y, dy)) 


def ln_Sinposterior(p, t, y, dy):
			ln_p = ln_Sinprior(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_Sinlikelihood(p, t, y, dy)
			return ln_l + ln_p



if (fmin_Fit):
### FIT SOURCE
	if (Fit_Src):
		Shell_File = "Source_Fit"
		param_names = ["Lfac", "beta", "phase", "incl"]
		No_Prd = True
		psrc = [5.79860006e-02,   6.83846646e-02,   3.36574394e-01, 9.66557496e-05]  # best fit for PG 1302
		fminFsrc_opt = sc.optimize.fmin(Fsrc_Err2,   psrc, args=(tsrt, Lumsrt, sigLsrt), full_output=1, disp=False,ftol=0.0001)[0]
		print fminFsrc_opt
		ShW1_p_opt = fminFsrc_opt
		ShW2_p_opt = fminFsrc_opt

	if (SinFit):
		Shell_File = "Sin_W1W2fmin_xefix"
		param_names = ["Amp", "phase", "Mag0"]
		print "Fmin optimizing W1"
		W1_sin_p_opt  = sc.optimize.fmin(SinErr2,     Sinp0_W1, args=(t_avg, W1_avg, W1_avsg), full_output=1, disp=False,ftol=0.0001)[0]
		print "Fmin optimizing W2"
		W2_sin_p_opt  = sc.optimize.fmin(SinErr2,     Sinp0_W2, args=(t_avg, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.0001)[0]
		ShW1_p_opt = W1_sin_p_opt
		ShW2_p_opt = W2_sin_p_opt


	if (ShellFit):
		Shell_File = "W1andW2fmin_Shell_xefix"
		param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$R_in$', r'$n_0$']
		print "Fmin optimizing W1"
		ShW1_p_opt  = sc.optimize.fmin(Shell_RegErr2,     ShW1_p0_0, args=(t_avg, W1args, RHS_table, T_table, W1_avg, W1_avsg, rem_is_Rin), full_output=1, disp=False,ftol=0.01)[0]
		print "Fmin optimizing W2"
		ShW2_p_opt  = sc.optimize.fmin(Shell_RegErr2,     ShW2_p0_0, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg, rem_is_Rin), full_output=1, disp=False,ftol=0.01)[0]
###BOTH
	if (fit_both):
		if (rem_is_Rin):
			ShW1_p_opt  = sc.optimize.fmin(ShellBoth_RegErr2,     Shboth_p0_0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg, rem_is_Rin), full_output=1, disp=False,ftol=0.01)[0]
			ShW2_p_opt  = ShW1_p_opt 
			Shell_File = "W1W2fmin_BOTH_SAMEreqRin_Shell_xefix"
			param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$Rin$', r'$n_0$']
			p1both = ShW1_p_opt
			p2both = ShW1_p_opt
		else:
			#ShW1_p_opt  = sc.optimize.fmin(ShellBoth_RegErr2,     ShW1_p0_0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg, rem_is_Rin), full_output=1, disp=False,ftol=0.01)[0]
			#ShW2_p_opt  = sc.optimize.fmin(ShellBoth_RegErr2,     ShW2_p0_0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg, rem_is_Rin), full_output=1, disp=False,ftol=0.01)[0]
			ShW1_p_opt  = sc.optimize.fmin(ShellBoth_RegErr2,     Shboth_p0_0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg, rem_is_Rin), full_output=1, disp=False,ftol=0.01)[0]
			Shell_File = "W1W2fmin_BOTH_rem_diff_Rin_Shell_xefix"
			param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$rem1$', r'$rem2$', r'$Rin$', r'$n_0$']
			p1both = [ShW1_p_opt[0], ShW1_p_opt[1], ShW1_p_opt[2], ShW1_p_opt[4], ShW1_p_opt[5]]
			p2both = [ShW1_p_opt[0], ShW1_p_opt[1], ShW1_p_opt[3], ShW1_p_opt[4], ShW1_p_opt[5]]

	if (ThickFit):
		if (W1fit):
			Shell_File = "W1_5p_fmin_Thick_Fsrcoffset_xefix"
			#p0 = [cosJ, costheta_T, Rin, p, n0]
			param_names = [r'sin($J$)',r'cos($\theta_T$)',r'$R_{in}$',r'$p$', r'$n_0$']
			print "Fmin optimizing W1"
			ShW1_p_opt  = sc.optimize.fmin(Thick_RegErr2_fmin,     ShW1_p0_0, args=(t_avg, W1args, RHS_table, T_table, W1_avg, W1_avsg), full_output=1, disp=False,ftol=0.1)[0]
			print "Fmin optimizing W2"
			ShW2_p_opt = ShW1_p_opt#ShW2_p_opt  = sc.optimize.fmin(Thick_RegErr2,     ShW2_p0_0, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.01)[0]
		if (W2fit):
			Shell_File = "W2_5p_fmin_ThickFsrcoffset_xefix"
			#p0 = [cosJ, costheta_T, Rin, p, n0]
			param_names = [r'sin($J$)',r'cos($\theta_T$)',r'$R_{in}$',r'$p$', r'$n_0$']
			print "Fmin optimizing W2"
			ShW2_p_opt  = sc.optimize.fmin(Thick_RegErr2_fmin,     ShW2_p0_0, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.1)[0]
			print "Fmin optimizing W1"
			ShW1_p_opt = ShW2_p_opt#ShW2_p_opt  = sc.optimize.fmin(Thick_RegErr2,     ShW2_p0_0, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.01)[0]
	if (NoFit):
		Shell_File = Shell_File + "_NoFit"
		if (rem_is_Rin):
			param_names = [r'cos($J$)',r'cos($\theta_T$)',r'$Rin$', r'$n_0$']
		else:	
			param_names = [r'cos($J$)',r'cos($\theta_T$)',r'$rem$',r'$Rin$', r'$n_0$']
		ShW1_p_opt = ShW1_p0_0
		ShW2_p_opt = ShW2_p0_0

	filename = "fmin"+Shell_File+"_results.txt"
	print "Printing Results"
	target = open(filename, 'w')
	#target.truncate()

		
	for i,name in enumerate(param_names):
		target.write("W1: {name}: {0:.4f}".format(ShW1_p_opt[i], name=name))
		target.write("\n")

	target.write("\n")

	for i,name in enumerate(param_names):
		target.write("W2: {name}: {0:.4f}".format(ShW2_p_opt[i], name=name))
		target.write("\n")
		
			

	target.close

	filename = "fmin"+Shell_File+"_results.txt"
	print "Printing Results"
	target = open(filename, 'w')
	#target.truncate()

		
	for i,name in enumerate(param_names):
		target.write("W1: {name}: {0:.4f}".format(ShW1_p_opt[i], name=name))
		target.write("\n")
	target.write("\n")

	for i,name in enumerate(param_names):
		target.write("W2: {name}: {0:.4f}".format(ShW2_p_opt[i], name=name))
		target.write("\n")

		
			

	target.close






if (emcee_Fit):
	print "Running MCMC (!)..."
	### Run MCMC

	if (sinFit_Src):
		ndim = 3
		nwalkers = ndim*12
		Shell_File = "sinFsrcFit_NOPrd"
		Fsrc_param_names = ['Amp','phase','mag0']

		Fsrc_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Sinposterior, threads=NThread,args=(tsrt, Lumsrt, sigLsrt))
		Fsrc_p0  = np.array([0.1269,   1403.3821,   14.8330])
		
		Fsrc_walker_p0 = np.random.normal(Fsrc_p0, np.abs(Fsrc_p0)*1E-4, size=(nwalkers, ndim))
		

					
		clen = 2048#4096*2
		Fsrc_pos,_,_ = Fsrc_sampler.run_mcmc(Fsrc_walker_p0 , clen)

		

		print "SAVING THE PICKLE mmmmm"
		with open("../emcee_data/Pickles/PG1302_Fsrc_%iwalkers.pickle" %clen, "w") as f1:
			pickle.dump((Fsrc_sampler.chain, Fsrc_sampler.lnprobability), f1)


				


		### OPEN OUTPUT DATA
		with open("../emcee_data/Pickles/PG1302_Fsrc_%iwalkers.pickle" %clen) as f1:
			Fsrc_chain,Fsrc_lnprobs = pickle.load(f1)



		Fsrc_flatchain   = np.vstack(Fsrc_chain[:,clen/2:])
		Fsrc_flatlnprobs = np.vstack(Fsrc_lnprobs[:,clen/2:])

					
				
		Fsrc_opt  = Fsrc_flatchain[Fsrc_flatlnprobs.argmax()]















	#sampler = emcee.EnsembleSampler(walkers, ndim, ln_posterior, args=(tsrt, W1args, RHStable, Ttable, W1_mag, W1_sig))
	if (SinFit):

		if (No_Prd):
			ndim = 3
			Shell_File = "SinFit_NOPrd"
			param_names = ['Amp','phase','mag0']	
		else:
			ndim = 4
			Shell_File = "SinFit_wPrd"
			param_names = ['Amp','Prd','phase','mag0']	
			
		nwalkers = ndim*12

		W1_sin_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Sinposterior, threads=NThread,args=(t_avg, W1_avg, W1_avsg))
		W2_sin_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Sinposterior, threads=NThread,args=(t_avg, W2_avg, W2_avsg))


		#p0 = np.array(p0)
		W1_sin_p0 = np.array(Sinp0_W1)
		W2_sin_p0 = np.array(Sinp0_W2)
		W1_sin_walker_p0 = np.random.normal(W1_sin_p0, np.abs(W1_sin_p0)*1E-4, size=(nwalkers, ndim))
		W2_sin_walker_p0 = np.random.normal(W2_sin_p0, np.abs(W2_sin_p0)*1E-4, size=(nwalkers, ndim))

					
		clen = 2048#4096*2
		W1_sin_pos,_,_ = W1_sin_sampler.run_mcmc(W1_sin_walker_p0 , clen)

		W2_sin_pos,_,_ = W2_sin_sampler.run_mcmc(W2_sin_walker_p0 , clen)

		print "SAVING THE PICKLE mmmmm"
		with open("../emcee_data/Pickles/PG1302_W1_sin_%iwalkers.pickle" %clen, "w") as f1:
			pickle.dump((W1_sin_sampler.chain, W1_sin_sampler.lnprobability), f1)

		with open("../emcee_data/Pickles/PG1302_W2_sin_%iwalkers.pickle" %clen, "w") as f1:
			pickle.dump((W2_sin_sampler.chain, W2_sin_sampler.lnprobability), f1)

				


		### OPEN OUTPUT DATA
		with open("../emcee_data/Pickles/PG1302_W1_sin_%iwalkers.pickle" %clen) as f1:
			W1_sin_chain,W1_sin_lnprobs = pickle.load(f1)

		with open("../emcee_data/Pickles/PG1302_W2_sin_%iwalkers.pickle" %clen) as f2:
			W2_sin_chain,W2_sin_lnprobs = pickle.load(f2)


		W1_sin_flatchain = np.vstack(W1_sin_chain[:,clen/2:])
		W1_sin_flatlnprobs = np.vstack(W1_sin_lnprobs[:,clen/2:])

		W2_sin_flatchain = np.vstack(W2_sin_chain[:,clen/2:])
		W2_sin_flatlnprobs = np.vstack(W2_sin_lnprobs[:,clen/2:])

					
				
		W1_sin_p_opt  = W1_sin_flatchain[W1_sin_flatlnprobs.argmax()]
		W2_sin_p_opt  = W2_sin_flatchain[W2_sin_flatlnprobs.argmax()]
















	if (ShellFit or ThickFit):
		
		
		# if (mpi_it):
		# 	import sys
		# 	from emcee.utils import MPIPool

		# 	pool = MPIPool()
		# 	if not pool.is_master():
		# 		pool.wait()
		# 		sys.exit(0)
		# 	ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Thposterior, pool=pool, args=(t_avg, W1args, RHS_table, T_table, W1_avg, W1_avsg))
		# 	pool.close()
		# else:
		if (ShellFit):
			ndim = 4
			nwalkers = ndim*12
			param_names = [r'sin($J$)',r'cos($\theta_T$)', r'$R_in$', r'$n_0$']
			if (W2fit):
				Shell_File = "W2_Shell_xefix"
				ShW1_p0 = np.array(ShW2_p0_0)
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Shposterior, threads=NThread, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg, rem_is_Rin))
			if (W1fit):
				Shell_File = "W1_Shell_xefix"
				ShW1_p0 = np.array(ShW1_p0_0)
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Shposterior, threads=NThread, args=(t_avg, W1args, RHS_table, T_table, W1_avg, W1_avsg, rem_is_Rin))
		if (fit_both):
			ShW1_p0 = np.array(Shboth_p0_0)
			if (rem_is_Rin):
				Shell_File = "_fitbothW1W2_rem_is_Rin"
				param_names = [r'sin($J$)',r'cos($\theta_T$)', r'$R_in$', r'$n_0$']
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_ShBothPosterior, threads=NThread, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg, rem_is_Rin))
			else:
				param_names = [r'cos($J$)',r'cos($\theta_T$)',r'$r_1$',r'$r_2$', r'$R_{in}$', r'$n_0$']
				ndim = 6
				nwalkers = ndim*8
				Shell_File = "_fitbothW1W2_rem1_rem2"
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_ShBothPosterior, threads=NThread, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg, rem_is_Rin))
			# else:
			# 	print "must choose W1 or W1 or both"
			# 	import sys
			# 	sys.exit(0)
		if (ThickFit):
			ndim = 4
			nwalkers = ndim*12
			param_names = [r'sin($J$)',r'cos($\theta_T$)',r'$R_{in}$', r'$p$', r'$n_0$']
			if (W2fit):
				Shell_File = "W2_Thick_xefix"
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Thposterior, threads=NThread, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg))
			elif (W1fit):
				Shell_File = "W1_Thick_xefix"
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Thposterior, threads=NThread, args=(t_avg, W1args, RHS_table, T_table, W1_avg, W1_avsg))
			else:
				print "must choose W1 or W1 to fit too (do both later)"
				import sys
				sys.exit(0)
		
		# if (W2fit):
		# 	ShW1_p0 = np.array(ShW2_p0_0)
		# elif (W1fit):
		# 	ShW1_p0 = np.array(ShW1_p0_0)
		# elif (fit_both):
		# 	ShW1_p0 = np.array(Shboth_p0_0)
		# else:
		# 	print "must choose W1 or W1 or both"
		# 	import sys
		# 	sys.exit(0)

		
		ShW1_walker_p0 = np.random.normal(ShW1_p0, np.abs(ShW1_p0)*1E-4, size=(nwalkers, ndim))
		
		#ShW2_p0 = np.array(ShW1_p0)
		#ShW2_walker_p0 = np.random.normal(ShW1_p0, np.abs(ShW1_p0)*1E-2, size=(nwalkers, ndim))
		

					
		clen = 1024
		ShW1_pos,_,_ = ShW1_sampler.run_mcmc(ShW1_walker_p0 , clen)


		print "SAVING THE PICKLE mmmmm"
		with open("../emcee_data/Pickles/PG1302_Shell_IRLE_"+Shell_File+"_%iwalkers.pickle" %clen, "w") as f1:
			pickle.dump((ShW1_sampler.chain, ShW1_sampler.lnprobability), f1)


		### OPEN OUTPUT DATA
		with open("../emcee_data/Pickles/PG1302_Shell_IRLE_"+Shell_File+"_%iwalkers.pickle" %clen) as f1:
			ShW1_chain,ShW1_lnprobs = pickle.load(f1)



		ShW1_flatchain = np.vstack(ShW1_chain[:,clen/2:])
		ShW1_flatlnprobs = np.vstack(ShW1_lnprobs[:,clen/2:])

					
				
		ShW1_p_opt  = ShW1_flatchain[ShW1_flatlnprobs.argmax()]
		if (fit_both):
			if (rem_is_Rin):
				p1both = ShW1_p_opt
				p2both = ShW1_p_opt
			else:
				p1both = [ShW1_p_opt[0], ShW1_p_opt[1], ShW1_p_opt[2], ShW1_p_opt[4], ShW1_p_opt[5]]
				p2both = [ShW1_p_opt[0], ShW1_p_opt[1], ShW1_p_opt[3], ShW1_p_opt[4], ShW1_p_opt[5]]





	print "ANALYSING MCMC (!)..."
	### MCMC ANALYSIS

	if (sinFit_Src):
		with open("../emcee_data/Pickles/PG1302_Fsrc_%iwalkers.pickle" %clen) as f1:
			Fsrc_chain,Fsrc_lnprobs = pickle.load(f1)
		


		##PLOT dem WALKERS
		for k in range(Fsrc_chain.shape[2]):
			plt.figure()
			#plt.figure()
			for i in range(Fsrc_chain.shape[0]):
				plt.plot(Fsrc_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
				plt.ylabel(Fsrc_param_names[k])
				plt.xlabel('steps')
			plt.savefig('../emcee_data/Fsrc_PG1302_%s_%iwalkers.png' %(Fsrc_param_names[k],clen))
			plt.clf()

		



		###CORNER PLOT	
		Fsrc_flatchain = np.vstack(Fsrc_chain[:,clen/2:])
		Fsrc_flatlnprobs = np.vstack(Fsrc_lnprobs[:,clen/2:])

				

		#import triangle
		import corner as triangle
		Fsrc_fig = triangle.corner(Fsrc_flatchain, labels=Fsrc_param_names)			
		Fsrc_fig.savefig('../emcee_data/Fsrc_PG1302_Corner_Plot_%iwalkers.png' %clen)




		## Do some stats on the walkers
		from scipy.stats import scoreatpercentile as scoretpercentile

		## max posterior + percentiles
		Fsrc_MAP_vals = Fsrc_flatchain[Fsrc_flatlnprobs.argmax()]
		Fsrc_perc = scoretpercentile(Fsrc_flatchain, [15,85], axis=0)


		filename = "Fsrc_results_"+Shell_File+"%iwalkers.txt" %clen
		print "Printing Results"
		target = open(filename, 'w')
		target.truncate()


		for i,name in enumerate(Fsrc_param_names):
			Fsrc_diff_minus = Fsrc_MAP_vals[i] - Fsrc_perc[0,i]
			Fsrc_diff_plus = Fsrc_perc[1,i] - Fsrc_MAP_vals[i]
			target.write("Fsrc: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(Fsrc_MAP_vals[i], Fsrc_diff_plus, Fsrc_diff_minus, name=name))
			target.write("\n")

		

		if (No_Prd):
			target.write("Period = %4g" %SinPrd)	
			target.write("\n")		
		
		
		
		Fsrc_mxprbs = zeros(nwalkers)
					

		for i in range(nwalkers):
			Fsrc_mxprbs[i] = max(Fsrc_lnprobs[i])
			
		
		chi2_pdf_Fsrc = -max(Fsrc_mxprbs)/(len(tsrt) - len(Fsrc_param_names) - 1)
		
		target.write("\n")		
		target.write("Fsrc_fit reduced chi2 =  %04g" %chi2_pdf_Fsrc)
		

			

		target.close












	if (SinFit):
					
		with open("../emcee_data/Pickles/PG1302_W1_sin_%iwalkers.pickle" %clen) as f1:
			W1_sin_chain,W1_sin_lnprobs = pickle.load(f1)
		with open("../emcee_data/Pickles/PG1302_W2_sin_%iwalkers.pickle" %clen) as f2:
			W2_sin_chain,W2_sin_lnprobs = pickle.load(f2)


		##PLOT dem WALKERS
		for k in range(W1_sin_chain.shape[2]):
			plt.figure()
			#plt.figure()
			for i in range(W1_sin_chain.shape[0]):
				plt.plot(W1_sin_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
				plt.ylabel(param_names[k])
				plt.xlabel('steps')
			plt.savefig('../emcee_data/W1_sin_PG1302_%s_%iwalkers.png' %(param_names[k],clen))
			plt.clf()

		for k in range(W2_sin_chain.shape[2]):
			plt.figure()
			#plt.figure()
			for i in range(W2_sin_chain.shape[0]):
				plt.plot(W2_sin_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
				plt.ylabel(param_names[k])
				plt.xlabel('steps')
			plt.savefig('../emcee_data/W2_sin_PG1302_%s_%iwalkers.png' %(param_names[k],clen))
			plt.clf()



		###CORNER PLOT	
		W1_sin_flatchain = np.vstack(W1_sin_chain[:,clen/2:])
		W1_sin_flatlnprobs = np.vstack(W1_sin_lnprobs[:,clen/2:])

		W2_sin_flatchain = np.vstack(W2_sin_chain[:,clen/2:])
		W2_sin_flatlnprobs = np.vstack(W2_sin_lnprobs[:,clen/2:])

				

		#import triangle
		import corner as triangle
		W1_sin_fig = triangle.corner(W1_sin_flatchain, labels=param_names)			
		W1_sin_fig.savefig('../emcee_data/W1_sin_PG1302_Corner_Plot_%iwalkers.png' %clen)


		W2_sin_fig = triangle.corner(W2_sin_flatchain, labels=param_names)			
		W2_sin_fig.savefig('../emcee_data/W2_sin_PG1302_Corner_Plot_%iwalkers.png' %clen)



		## Do some stats on the walkers
		from scipy.stats import scoreatpercentile as scoretpercentile

		# ## mean +- stdev
		# #Bmean_vals = np.mean(flatchain, axis=0)
		# #Bstddev = np.std(flatchain, axis=0)


		# #for i,name in enumerate(param_names):
		# #	print("{name}: {0:.3f} +/- {1:3f}".format(mean_vals[i], stddev[i], name=name))

		# ## median + percentiles
		# med_vals = np.median(flatchain, axis=0)
		# perc = scoreatpercentile(flatchain, [15,85], axis=0)

		# for i,name in enumerate(param_names):
		# 	diff_minus = med_vals[i] - perc[0,i]
		# 	diff_plus = perc[1,i] - med_vals[i]
		# 	print("{name}: {0:.3f} + {1:.3f} - {2:.3f}".format(med_vals[i], diff_plus, diff_minus, name=name))


		## max posterior + percentiles
		W1_sin_MAP_vals = W1_sin_flatchain[W1_sin_flatlnprobs.argmax()]
		W1_sin_perc = scoretpercentile(W1_sin_flatchain, [15,85], axis=0)

		W2_sin_MAP_vals = W2_sin_flatchain[W2_sin_flatlnprobs.argmax()]
		W2_sin_perc = scoretpercentile(W2_sin_flatchain, [15,85], axis=0)




		# for i,name in enumerate(param_names):
		# 	W1_sin_diff_minus = W1_sin_MAP_vals[i] - W1_sin_perc[0,i]
		# 	W1_sin_diff_plus = W1_sin_perc[1,i] - W1_sin_MAP_vals[i]
		# 	print("W1: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(W1_sin_MAP_vals[i], W1_sin_diff_plus, W1_sin_diff_minus, name=name))

		# for i,name in enumerate(param_names):
		# 	W2_sin_diff_minus = W2_sin_MAP_vals[i] - W2_sin_perc[0,i]
		# 	W2_sin_diff_plus = W2_sin_perc[1,i] - W2_sin_MAP_vals[i]
		# 	print("W2: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(W2_sin_MAP_vals[i], W2_sin_diff_plus, W2_sin_diff_minus, name=name))
					
		# W1_sin_mxprbs = zeros(nwalkers)
		# W2_sin_mxprbs = zeros(nwalkers)
					
		# for i in range(nwalkers):
		# 	W1_sin_mxprbs[i] = max(W1_sin_lnprobs[i])
		# 	W2_sin_mxprbs[i] = max(W2_sin_lnprobs[i])

					
		# print "W1 Max LnP = ", max(W1_sin_mxprbs)
		# print "W2 Max LnP = ", max(W2_sin_mxprbs)		

		filename = "Sin_results_"+Shell_File+"%iwalkers.txt" %clen
		print "Printing Results"
		target = open(filename, 'w')
		target.truncate()


		for i,name in enumerate(param_names):
			W1_sin_diff_minus = W1_sin_MAP_vals[i] - W1_sin_perc[0,i]
			W1_sin_diff_plus = W1_sin_perc[1,i] - W1_sin_MAP_vals[i]
			target.write("W1: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(W1_sin_MAP_vals[i], W1_sin_diff_plus, W1_sin_diff_minus, name=name))
			target.write("\n")

		for i,name in enumerate(param_names):
			W2_sin_diff_minus = W2_sin_MAP_vals[i] - W2_sin_perc[0,i]
			W2_sin_diff_plus = W2_sin_perc[1,i] - W2_sin_MAP_vals[i]
			target.write("W2: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(W2_sin_MAP_vals[i], W2_sin_diff_plus, W2_sin_diff_minus, name=name))
			target.write("\n")	

		if (No_Prd):
			target.write("Period = %4g" %SinPrd)	
			target.write("\n")		
		
		
		
		W1_sin_mxprbs = zeros(nwalkers)
		W2_sin_mxprbs = zeros(nwalkers)
					

		for i in range(nwalkers):
			W1_sin_mxprbs[i] = max(W1_sin_lnprobs[i])
			W2_sin_mxprbs[i] = max(W2_sin_lnprobs[i])
		
		chi2_pdf_W1 = -max(W1_sin_mxprbs)/(len(W1_avg) - len(param_names) - 1)
		chi2_pdf_W2 = -max(W2_sin_mxprbs)/(len(W2_avg) - len(param_names) - 1)
		target.write("\n")		
		target.write("Shell W1 reduced chi2 =  %04g" %chi2_pdf_W1)
		target.write("\n")
		target.write("Shell W2 reduced chi2 =  %04g" %chi2_pdf_W2)

			

		target.close





	if (ShellFit or ThickFit):
		#param_names = ['beta','cosJ','Rin','thetT', 'n0']			
			
		


		with open("../emcee_data/Pickles/PG1302_Shell_IRLE_"+Shell_File+"_%iwalkers.pickle" %clen) as f1:
			ShW1_chain,ShW1_lnprobs = pickle.load(f1)
		


		##PLOT dem WALKERS
		for k in range(ShW1_chain.shape[2]):
			plt.figure(param_names[k])
			#plt.figure()
			for i in range(ShW1_chain.shape[0]):
				plt.plot(ShW1_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
				plt.ylabel(param_names[k])
				plt.xlabel('steps')
			plt.savefig("../emcee_data/"+Shell_File+"_PG1302_%s_%iwalkers.png" %(param_names[k],clen))


		



		###CORNER PLOT	
		ShW1_flatchain = np.vstack(ShW1_chain[:,clen/2:])
		ShW1_flatlnprobs = np.vstack(ShW1_lnprobs[:,clen/2:])
			

		#import triangle
		import corner as triangle
		ShW1_fig = triangle.corner(ShW1_flatchain, labels=param_names)			
		ShW1_fig.savefig("../emcee_data/"+Shell_File+"_PG1302_Corner_Plot_%iwalkers.png" %clen)


		## Do some stats on the walkers
		from scipy.stats import scoreatpercentile as scoretpercentile

		
		## max posterior + percentiles
		ShW1_MAP_vals = ShW1_flatchain[ShW1_flatlnprobs.argmax()]
		ShW1_perc = scoretpercentile(ShW1_flatchain, [15,85], axis=0)


		filename = "SHW1_results_"+Shell_File+"%iwalkers.txt" %clen
		print "Printing Results"
		target = open(filename, 'w')
		target.truncate()

		
		for i,name in enumerate(param_names):
			ShW1_diff_minus = ShW1_MAP_vals[i] - ShW1_perc[0,i]
			ShW1_diff_plus = ShW1_perc[1,i] - ShW1_MAP_vals[i]
			target.write("W1: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(ShW1_MAP_vals[i], ShW1_diff_plus, ShW1_diff_minus, name=name))
			target.write("\n")
		
		ShW1_mxprbs = zeros(nwalkers)
					
		for i in range(nwalkers):
			ShW1_mxprbs[i] = max(ShW1_lnprobs[i])
		

		target.write("\n")		
		target.write("Shell W1 Max LnP =  %04g" %max(ShW1_mxprbs))

			

		target.close






#####------PLOT SOLUTION------####
### PLOT POINTS


#tsrt = tsrt #- 49100
#t_MJD = t_MJD #- 49100

Nt=40

ttopt = np.linspace(tsrt[0]-100, t_MJD[len(t_MJD)-1]+100,       Nt)


#opti = -2.5*np.log10(Fsrc(ttopt*3600.*24, Dst, ma.pi/2., 0.0, Lav, betst, Inc, OmPG, 1.1)/FVbndRel)
if (Fit_Src):
	#opti = -2.5*np.log10(Fsrc((ttopt*3600.*24 - fminFsrc_opt[2]*2.*ma.pi/OmPG), Dst, ma.pi/2., 0.0, fminFsrc_opt[0]*Lav, fminFsrc_opt[1], np.arccos(0.067/fminFsrc_opt[1]), OmPG, 1.1)/FVbndRel)
	opti = -2.5*np.log10(Fsrc((ttopt*3600.*24 - 0.0), Dst, ma.pi/2., 0.0, fminFsrc_opt[0]*Lav, fminFsrc_opt[1], fminFsrc_opt[3], OmPG, 1.1)/FVbndRel)
#elif (sinFit_Src):
#	opti =  plt.plot(ttopt, sinPoint(Fsrc_opt, (ttopt+50000)/(1.+zPG1302)), linestyle = '--', color='orange', linewidth=2)
else:
	fminFsrc_opt = [ 5.79860006e-02,   6.83846646e-02,   3.36574394e-01, 9.66557496e-05]
	#opti = -2.5*np.log10(Fsrc((ttopt*3600.*24 - fminFsrc_opt[2]*2.*ma.pi/OmPG), Dst, ma.pi/2., 0.0, fminFsrc_opt[0]*Lav, fminFsrc_opt[1], np.arccos(0.067/fminFsrc_opt[1]), OmPG, 1.1)/FVbndRel)
	opti = -2.5*np.log10(Fsrc((ttopt*3600.*24 - 0.0), Dst, ma.pi/2., 0.0, fminFsrc_opt[0]*Lav, fminFsrc_opt[1], fminFsrc_opt[3], OmPG, 1.1)/FVbndRel)


ttopt = (ttopt*(1.+zPG1302) - 50000)
t_avg = (t_avg*(1.+zPG1302) - 50000)
tsrt  = (tsrt*(1.+zPG1302) - 50000)
t_MJD = (t_MJD*(1.+zPG1302) - 50000)

## z adjust time parameers
## PERIOD
#W1_sin_p_opt[1] = W1_sin_p_opt[1]*(1.+zPG1302) 
#W2_sin_p_opt[1] = W2_sin_p_opt[1]*(1.+zPG1302) 

print "PLOTTING"

plt.figure()
plt.errorbar(tsrt, Lumsrt-3.0, yerr=sigL, linestyle="none", color = "blue", alpha=0.5) #alpha=0.1
if (sinFit_Src):
	Fs =  plt.plot(ttopt, sinPoint(Fsrc_opt, (ttopt+50000)/(1.+zPG1302))-3.0, linestyle = '--', color='blue', linewidth=2)
else:
	Fs = plt.plot(ttopt, opti-3.0, linestyle = '--', color='blue', linewidth=2)


W1dat   = plt.errorbar(t_MJD, W1_mag, yerr=W1_sig, linestyle="none", color='orange', alpha=1., elinewidth=1.5)
W2dat   = plt.errorbar(t_MJD, W2_mag+0.5, yerr=W2_sig, linestyle="none", color='red', alpha=1., elinewidth=1.5)

W1av   = plt.errorbar(t_avg, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)
W2av   = plt.errorbar(t_avg, W2_avg+0.5, yerr=W2_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)


if (SinFit):
	W1sinsoln = plt.plot(ttopt, sinPoint(W1_sin_p_opt, (ttopt+50000)/(1.+zPG1302)), linestyle = '--', color='orange', linewidth=2)
	W2sinsoln = plt.plot(ttopt, sinPoint(W2_sin_p_opt, (ttopt+50000)/(1.+zPG1302))+0.5, linestyle = '--', color='red', linewidth=2)
if (ShellFit):
	if (fmin_Fit):
		W1shell = plt.plot(ttopt, magPoint_Shell(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table, rem_is_Rin), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_Shell(ShW2_p_opt, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table, rem_is_Rin)+0.5, linestyle = '--', color='red', linewidth=2)
	else:
		W1shell = plt.plot(ttopt, magPoint_Shell(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table, rem_is_Rin), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_Shell(ShW2_p0_0, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table, rem_is_Rin)+0.5, linestyle = '--', color='red', linewidth=2)
if (fit_both or pltboth):
		W1shell = plt.plot(ttopt, magPoint_Shell(p1both, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table, rem_is_Rin), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_Shell(p2both, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table, rem_is_Rin)+0.5, linestyle = '--', color='red', linewidth=2)
if (ThickFit):	
	if (fmin_Fit):
		W1shell = plt.plot(ttopt, magPoint_Thick(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_Thick(ShW2_p_opt, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
	else:
		W1shell = plt.plot(ttopt, magPoint_Thick(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_Thick(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
if (NoFit and pltShell):
	W1shell = plt.plot(ttopt, magPoint_Shell(ShW1_p0_0, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table, rem_is_Rin), linestyle = '--', color='orange', linewidth=2)
	W2shell = plt.plot(ttopt, magPoint_Shell(ShW1_p0_0, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table, rem_is_Rin)+0.5, linestyle = '--', color='red', linewidth=2)
if (NoFit and pltThick):
	W1shell = plt.plot(ttopt, magPoint_Thick(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
	W2shell = plt.plot(ttopt, magPoint_Thick(ShW2_p_opt, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
	


plt.grid(b=True, which='both')
		#plt.legend( [ s1[0], IR1[0], s2[0], IR2[0], s3[0], IR3[0]  ], (r'$i=0$','',   r'$i=\pi/4$','',   r'$i=\pi/2$', ''), loc='upper right')

plt.xlabel(r"$t$ [MJD]")
plt.ylabel("mag")
#plt.xlim(52000, 57500)
plt.xlim(3000, max(ttopt))
#plt.ylim(10.5, 11.5)
#plt.ylim(plt.ylim(10.5, 12.3)[::-1])

		#plt.show()
plt.savefig("../emcee_data/"+Shell_File+"BestFit.png")
plt.clf()

print "ALL DONE"


