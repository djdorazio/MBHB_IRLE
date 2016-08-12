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
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams.update({'font.size': 16})
import matplotlib.pyplot as plt

# from scipy import optimize
# from scipy.optimize import fmin

from FluxFuncs_IRLE import *
from ErrFuncs_IRLE import *
from emcee_Funcs import *

################################
###############################
### OPTIONS
################################
################################
Shell_OptThin = True
Fit_mag0 = True
TwoRs = False
fmin_start = False

##multiprocessing
NThread = 32

#Temp table resolution
NTemp = 1800
Tmin = 100.
Tsub  = 1800.


################################
###############################
### Define Constants
################################
################################
nne = 0.0
nu0 = numicron#/1.5

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
aeff = (c/nu0)/(2.*ma.pi) #(1 micrometer wavelength /(2pi) )


Dst = 1.4*10**9*pc2cm
Rrout = 1.0*Rde # DEFUNCT PARAMETER



## Wise band numbers
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3



nuVbnd = c/(5.45*10**(-5))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-21)*8.8560*10**(13)#(W1mn + W1mx)/2
FW2Rel = 1.71787*10**(-21)*6.4451*10**(13)#(W2mn + W2mx)/2










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
##TABULATE T's and RHSs
print "Creating Temp look up tables..."
RHS_table = np.zeros(NTemp)
T_table = np.linspace(Tmin, Tsub, NTemp)
for i in range(NTemp):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)


















if (Shell_OptThin):
	if (fmin_start):
		print "Running fmin to start"
		from scipy import optimize
		from scipy.optimize import fmin

		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, pp, Rrout, nu0, nne, betst] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, pp, Rrout, nu0, nne, betst] 

		if (TwoRs):
			p0 = [ 0.94640088,  0.12592987,  2.96947227,  4.41782365]
			popt  = sc.optimize.fmin(DOP_OpThin_TorShell_Err2_TwoRs,    p0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False, ftol=0.01)[0]
	
		if (Fit_mag0):
			p0 = [ 0.0121,  0.1168,  1.2238,  0.5365]
			#0.97892352  0.12836802  4.13484639 -0.20489902 chi 254
			popt  = sc.optimize.fmin(OpThin_TorShell_Err2_mag0,    p0, args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg), full_output=1, disp=False, ftol=0.01)[0]
	


	print "SETTING UP OPT-THIN GEO-THIN TORUS SHELL MCMC (!)..."
	if (TwoRs):
		Shell_File = "DOP_GeoThin_OptThin_NTemp%g_Tsub%g_" %(NTemp, Tsub)
		param_names = [r'$\sin{J}$', r'$\cos{\theta_T}$', r'$R_1$', r'$R_2$']
		###MEASURED VALUES (edge on ring):
		p0 = [0.01, 0.125, 1.75, 1.75]
	else:
		if (Fit_mag0):
			Shell_File = "DOP_GeoThin_OptThin_Fitmag0W1_NTemp%g_Tsub%g_TwoRs_" %(NTemp, Tsub)
			param_names = [r'$\sin{J}$', r'$\cos{\theta_T}$', r'$R_{\rm{d}}$', r'$mag^{\rm{W1}}_0$']
			###MEASURED VALUES (edge on ring):	
			p0 = [0.0121,  0.1168,  1.2238,  0.5365]
		else:
			Shell_File = "DOP_GeoThin_OptThin_NTemp%g_Tsub%g_TwoRs_" %(NTemp, Tsub)
			param_names = [r'$\sin{J}$', r'$\cos{\theta_T}$', r'$R_{\rm{d}}$']
			###MEASURED VALUES (edge on ring):	
			p0 = [0.01, 0.125, 1.75]

	ndim = len(param_names)	
	nwalkers = ndim*8

	#Best fit from fmin (Dop_Fitting.py)
	#p0 = [-0.85664946,  0.87693384,  9.48441709, -1.3605258]
	#p0 = [0.7,  0.7,  0.5]
	#Best for from 2048 36 MCMC chi2 ~100, but too dim when integrate at better resolution - better resolution aslo gets rid of large amplitude as expected
	#p0 = [0.2065, 0.0735, 10.3518]
	# make brighter
	#p0 = [0.2065, 0.2, 10.3518]


	###MEASURED VALUES (edge on ring):
	#
	p0 = np.array(p0)
	## args of non changing parameters to pass 
	##Rrout is no longer used
	W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, pp, Rrout, nu0, nne, betst] 
	W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, pp, Rrout, nu0, nne, betst] 

	if (TwoRs):
		sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_SHThin_posterior_TwoRs, threads=NThread,args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg))
	else:
		if (Fit_mag0):
			sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_SHThin_posterior_mag0, threads=NThread,args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg))
		else:
			sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_SHThin_posterior, threads=NThread,args=(t_avg, W1args, W2args, RHS_table, T_table, W1_avg, W1_avsg, W2_avg, W2_avsg))


	walker_p0 = np.random.normal(p0, np.abs(p0)*1E-4, size=(nwalkers, ndim))


	clen = 2048
	pos,_,_ = sampler.run_mcmc(walker_p0 , clen)


	print "SAVING THE PICKLE mmmmm"
	with open("../emcee_data/Pickles/PG1302_Shell_IRLE_"+Shell_File+"_%iwalkers.pickle" %clen, "w") as f1:
		pickle.dump((sampler.chain, sampler.lnprobability), f1)


	### OPEN OUTPUT DATA
	with open("../emcee_data/Pickles/PG1302_Shell_IRLE_"+Shell_File+"_%iwalkers.pickle" %clen) as f1:
		chain,lnprobs = pickle.load(f1)



	flatchain   = np.vstack(chain[:,clen/2:])
	flatlnprobs = np.vstack(lnprobs[:,clen/2:])


	p_opt  = flatchain[flatlnprobs.argmax()]




	print "ANALYSING MCMC (!)..."
	### MCMC ANALYSIS

	with open("../emcee_data/Pickles/PG1302_Shell_IRLE_"+Shell_File+"_%iwalkers.pickle" %clen) as f1:
		chain,lnprobs = pickle.load(f1)
	


	##PLOT dem WALKERS
	for k in range(chain.shape[2]):
		plt.figure(param_names[k])
		#plt.figure()
		for i in range(chain.shape[0]):
			plt.plot(chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
			plt.ylabel(param_names[k])
			plt.xlabel('steps')
		plt.savefig("../emcee_data/"+Shell_File+"_PG1302_%s_%iwalkers.png" %(param_names[k],clen))

	###CORNER PLOT	
	flatchain   = np.vstack(chain[:,clen/2:])
	flatlnprobs = np.vstack(lnprobs[:,clen/2:])
	



	




	#import triangle
	import corner as triangle
	fig = triangle.corner(flatchain, labels=param_names)			
	fig.savefig("../emcee_data/"+Shell_File+"_PG1302_Corner_Plot_%iwalkers.png" %clen)


	## Do some stats on the walkers
	from scipy.stats import scoreatpercentile as scoretpercentile

	
	## max posterior + percentiles
	MAP_vals = flatchain[flatlnprobs.argmax()]
	perc = scoretpercentile(flatchain, [15,85], axis=0)


	filename = "MCMCresults_"+Shell_File+"%iwalkers.txt" %clen
	print "Printing Results"
	target = open(filename, 'w')
	target.truncate()

	
	for i,name in enumerate(param_names):
		diff_minus = MAP_vals[i] - perc[0,i]
		diff_plus = perc[1,i] - MAP_vals[i]
		target.write("W1: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(MAP_vals[i], diff_plus, diff_minus, name=name))
		target.write("\n")
	
	mxprbs = zeros(nwalkers)
				
	for i in range(nwalkers):
		mxprbs[i] = max(lnprobs[i])
	

	target.write("\n")		
	target.write("Dop Shell Max LnP =  %04g" %max(mxprbs))

		

	target.close()


	print "PLOTTING BEST FIT LIGHT CURVES"

	from Gen_Plot import *
	Plot_Shell_Thin_Dop(p_opt, Fit_mag0, TwoRs, 60,  Shell_File,  W1args, W2args, RHS_table, T_table,     tsrt, t_avg, t_MJD,    Lumsrt, W1_mag, W2_mag, W1_avg, W2_avg,   sigL, W1_sig, W2_sig,W1_avsg, W2_avsg)

else:


	W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, pp, Rrout, nu0, nne, betst] 
	W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, pp, Rrout, nu0, nne, betst] 

	if (TwoRs):
		Shell_File = "DOP_GeoThin_OptThin_NTemp%g_Tsub%g_" %(NTemp, Tsub)
		p0 = [0.01, 0.125, 1.75, 1.75]
	else:
		if (Fit_mag0):
			Shell_File = "DOP_GeoThin_OptThin_Fitmag0W1_NTemp%g_Tsub%g_" %(NTemp, Tsub)
			p0 = [0.0121,  0.1168,  1.2238,  0.5365]
		else:
			Shell_File = "DOP_GeoThin_OptThin_NTemp%g_Tsub%g_" %(NTemp, Tsub)
			p0 = [0.01, 0.125, 1.75]



	print "PLOTTING TEST LIGHT CURVES - DOPPLER"
	from Gen_Plot import *
	Plot_Shell_Thin_Dop(p0, Fit_mag0, TwoRs, 40,  Shell_File,  W1args, W2args, RHS_table, T_table,     tsrt, t_avg, t_MJD,    Lumsrt, W1_mag, W2_mag, W1_avg, W2_avg,   sigL, W1_sig, W2_sig,W1_avsg, W2_avsg)


















