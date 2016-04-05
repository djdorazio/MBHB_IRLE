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
NoFit = True
pltShell = False
pltThick = True

emcee_Fit = False
W1fit = False
W2fit = True

fmin_Fit = False
SinFit = False
ShellFit = True
ThickFit = False
## multiprocessing
NThread = 48
mpi_it = False
if (NoFit):
	emcee_Fit = False
	fmin_Fit = True

	ShellFit = False
	ThickFit = False
	SinFit   = False


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
OmPG = Omb*2.*np.pi/4.1
#alphnu = 0.0

Rorb = c*2.*np.pi/Omb
Ompc = 2.*np.pi*c/pc2cm/2.


## TEST VALUES
Lav = L0
betst = 0.08
Inc = ma.acos(0.067/betst)#0.*np.pi/4.
Ombn = OmPG
alph = -1.0

Rde = RdPG
pp = 2.0
thetTst = 1.*np.pi/4
JJt =4.*np.pi/8
aeff = 0.16*10**(-4) #(0.1 micrometer is an average ISM dust grain size)


Dst = 1.4*10**9*pc2cm
Rrout = 20.0*Rde

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


#p0 = [beta0, cosJJ0, Rin0, nDust0]
Sinp0_W1 = [0.1, 365*4.1, 1.0, 11.3]
Sinp0_W2 = [0.1, 365*4.1, 1.0, 10.3]
#[beta, cosJ, Rde(units of RdPG), theta_T, ndust(units of nDust0)]
#ShW1_p0  = [0.2,  1.,  1.4,  1.18427411,  1.7]
#[cosJ, Rde(units of RdPG) ndust(units of nDust0)]

#ShW1_p0_0  = [ 0.99455796,  1.6,  0.62607537]
#ShW2_p0_0  = [0.99935378,  2.0,  0.62607537]#0.81884671]


if (ShellFit):
	#p0 = [cosJ, costheta_T, Rin, n0]
	#ShW1_p0_0  = [ 0.001,  0.6905, 1.4392,  0.5880]
	#ShW2_p0_0  = [ 0.0009,  0.6035, 1.09,  2.5117]
	ShW1_p0_0  = [ 0.001,   0.6905, 1.4392,  0.5880] ## this fit give p~4
	ShW2_p0_0  = [ 0.0009,  0.6035, 1.0947,  2.5117]
	W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
	W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
if (ThickFit):
	#p0 = [cosJ, costheta_T, p, n0]
	ShW1_p0_0  = [ 0.0016,  0.9115, 0.5166,  0.1132]
	ShW2_p0_0  = [ 0.0016,  0.9115, 0.5166,  0.1132]
	W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, Rde, Rrout,  aeff, nu0, nne, betst] 
	W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, Rde, Rrout,  aeff, nu0, nne, betst] 
if (NoFit):
	#ShW1_p0_0  = [ 0.001,   0.6905, 1.4392,  0.5880]
	#ShW2_p0_0  = [ 0.0009,  0.6035, 1.0947,  2.5117]
	ShW1_p0_0  = [ 0.0016,  0.9115, 0.5166,  0.1132]
	ShW2_p0_0  = [ 0.0016,  0.9115, 0.5166,  0.1132]
	if (pltShell or pltThick):
		W1args = [FW1Rel, W1mn, W1mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 
		W2args = [FW2Rel, W2mn, W2mx, Dst, Lav, Ombn, alph, pp, Rrout,  aeff, nu0, nne, betst] 



#Targs = [Lav, betst, Inc, Ombn, alph, n0, Rde, pp, thetTst, JJt, aeff, nu0, nne]






print "Importing Data to fit..."
#Import Data to fit

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



t_avg = np.array(t_avg)
W1_avg = np.array(W1_avg)
W2_avg = np.array(W2_avg)
W1_avsg = np.array(W1_avsg)
W2_avsg = np.array(W2_avsg)

#import sys
#sys.exit(0)
### averaging data ^####

if (ShellFit or ThickFit or pltShell or pltThick):
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
def Shell_RegErr2(p, t, THEargs, RHStable, Ttable, y, dy):
	print "EVAL", p
	t1=time.clock()
	chi = (y - magPoint_Shell(p, t, THEargs, RHStable, Ttable)) / dy
	#nLnP = sum(chi*chi)
	t2=time.clock()
	#print(chi2)
	print(t2-t1)
	return sum(chi*chi)


def Thick_RegErr2(p, t, THEargs, RHStable, Ttable, y, dy):
	print "EVAL", p
	t1=time.clock()
	chi = (y - magPoint_Thick(p, t, THEargs, RHStable, Ttable)) / dy
	#nLnP = sum(chi*chi)
	t2=time.clock()
	#print(chi2)
	print(t2-t1)
	return sum(chi*chi)


def sinPoint(params, t):
	Amp, Prd, phs, mag0 = params
	#Amp, phs, mag0 = params
	#Prd=1884.
	return Amp*np.sin( (2.*ma.pi)/Prd*(t - Prd*phs) ) + mag0 


def SinErr2(p, t, y, dy):
	print "EVAL", p
	chi = (y - sinPoint(p, t) )/ dy
	return sum(chi*chi)
	#print(chi2)
	#return nLnP



### MCMC - Set up priors
def ln_prior(params):
			#beta, cosJJ, Rin, thetT, n0 = p
			cosJJ, cosTT, pp, n0 = params
			#if beta < 0.07 or beta > 0.5:
			#	return -np.inf
					
			if cosJJ < 0 or cosJJ > 1:
				return -np.inf

			if cosTT < 0 or cosTT > 1:
				return -np.inf
					
			if pp <= 1.0:
				return -np.inf
					

			if n0 <= 0.0:
				return -np.inf
					
				
			return 0.

def ln_Sinprior(p):
			Amp, Prd, phs, mag0 = p
			if Amp < 0.:
				return -np.inf
					
			if Prd < 0.:
				return -np.inf

			if phs < -1.0 or phs > 1.0:
				return -np.inf
					
				
			return 0.

### MCMC - Set up posteriors
def ln_Shlikelihood(p, t, Wargs, RHStable, Ttable, y, dy):
			return -(Shell_RegErr2(p, t, Wargs, RHStable, Ttable, y, dy)) #+ RegErr2(p, t, W2args, RHStable, Ttable, y2, dy2))

def ln_Thlikelihood(p, t, Wargs, RHStable, Ttable, y, dy):
			return -(Thick_RegErr2(p, t, Wargs, RHStable, Ttable, y, dy)) #+ RegErr2(p, t, W2args, RHStable, Ttable, y2, dy2))



def ln_Shposterior(p, t, Wargs, RHStable, Ttable, y, dy):
			ln_p = ln_prior(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_Shlikelihood(p, t, Wargs, RHStable, Ttable, y, dy)
			return ln_l + ln_p


def ln_Thposterior(p, t, Wargs, RHStable, Ttable, y, dy):
			ln_p = ln_prior(p)
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
	if (SinFit):
		print "Fmin optimizing W1"
		W1_sin_p_opt  = sc.optimize.fmin(SinErr2,     Sinp0_W1, args=(t_avg, W1_avg, W1_avsg), full_output=1, disp=False,ftol=0.0001)[0]
		print "Fmin optimizing W2"
		W2_sin_p_opt  = sc.optimize.fmin(SinErr2,     Sinp0_W2, args=(t_avg, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.0001)[0]

	if (ShellFit):
		Shell_File = "W1W2fmin_Shell"
		param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$R_in$', r'$n_0$']
		print "Fmin optimizing W1"
		ShW1_p_opt  = sc.optimize.fmin(Shell_RegErr2,     ShW1_p0_0, args=(t_avg, W1args, RHS_table, T_table, W1_avg, W1_avsg), full_output=1, disp=False,ftol=0.01)[0]
		print "Fmin optimizing W2"
		ShW2_p_opt  = sc.optimize.fmin(Shell_RegErr2,     ShW2_p0_0, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.01)[0]
	if (ThickFit):
		Shell_File = "W1fmin_Thick"
		param_names = [r'cos($J$)',r'cos($\theta_T$)',r'$p$', r'$n_0$']
		print "Fmin optimizing W1"
		ShW1_p_opt  = sc.optimize.fmin(Thick_RegErr2,     ShW1_p0_0, args=(t_avg, W1args, RHS_table, T_table, W1_avg, W1_avsg), full_output=1, disp=False,ftol=0.01)[0]
		print "Fmin optimizing W2"
		ShW2_p_opt = ShW1_p_opt#ShW2_p_opt  = sc.optimize.fmin(Thick_RegErr2,     ShW2_p0_0, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg), full_output=1, disp=False,ftol=0.01)[0]
	if (NoFit):
		Shell_File = "NoFit"
		param_names = [r'cos($J$)',r'cos($\theta_T$)',r'$p$', r'$n_0$']
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

				
	#sampler = emcee.EnsembleSampler(walkers, ndim, ln_posterior, args=(tsrt, W1args, RHStable, Ttable, W1_mag, W1_sig))
	if (SinFit):
		ndim = 4
		nwalkers = ndim*32

		W1_sin_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Sinposterior, threads=NThread,args=(t_avg, W1_avg, W1_avsg))
		W2_sin_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Sinposterior, threads=NThread,args=(t_avg, W2_avg, W2_avsg))


		#p0 = np.array(p0)
		W1_sin_p0 = np.array(Sinp0_W1)
		W2_sin_p0 = np.array(Sinp0_W2)
		W1_sin_walker_p0 = np.random.normal(W1_sin_p0, np.abs(W1_sin_p0)*1E-4, size=(nwalkers, ndim))
		W2_sin_walker_p0 = np.random.normal(W2_sin_p0, np.abs(W2_sin_p0)*1E-4, size=(nwalkers, ndim))

					
		clen = 512#4096*2
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
			param_names = [r'cos($J$)',r'cos($\theta_T$)', r'$R_in$', r'$n_0$']
			if (W2fit):
				Shell_File = "W2_Shell"
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Shposterior, threads=NThread, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg))
			elif (W1fit):
				Shell_File = "W1_Shell"
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Shposterior, threads=NThread, args=(t_avg, W1args, RHS_table, T_table, W1_avg, W1_avsg))
			else:
				print "must choose W1 or W1 to fit too (do both later)"
				import sys
				sys.exit(0)
		if (ThickFit):
			ndim = 4
			nwalkers = ndim*12
			param_names = [r'cos($J$)',r'cos($\theta_T$)',r'$p$', r'$n_0$']
			if (W2fit):
				Shell_File = "W2_Thick"
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Thposterior, threads=NThread, args=(t_avg, W2args, RHS_table, T_table, W2_avg, W2_avsg))
			elif (W1fit):
				Shell_File = "W1_Thick"
				ShW1_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Thposterior, threads=NThread, args=(t_avg, W1args, RHS_table, T_table, W1_avg, W1_avsg))
			else:
				print "must choose W1 or W1 to fit too (do both later)"
				import sys
				sys.exit(0)
		
		if (W2fit):
			ShW1_p0 = np.array(ShW2_p0_0)
		elif (W1fit):
			ShW1_p0 = np.array(ShW1_p0_0)
		else:
			print "must choose W1 or W1 to fit too (do both later)"
			import sys
			sys.exit(0)

		
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




	print "ANALYSING MCMC (!)..."
	### MCMC ANALYSIS

	#param_names = ['beta','cosJ','Rin','n0']

	if (SinFit):
		param_names = ['Amp','Prd','phase','mag0']	
					


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




		for i,name in enumerate(param_names):
			W1_sin_diff_minus = W1_sin_MAP_vals[i] - W1_sin_perc[0,i]
			W1_sin_diff_plus = W1_sin_perc[1,i] - W1_sin_MAP_vals[i]
			print("W1: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(W1_sin_MAP_vals[i], W1_sin_diff_plus, W1_sin_diff_minus, name=name))

		for i,name in enumerate(param_names):
			W2_sin_diff_minus = W2_sin_MAP_vals[i] - W2_sin_perc[0,i]
			W2_sin_diff_plus = W2_sin_perc[1,i] - W2_sin_MAP_vals[i]
			print("W2: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(W2_sin_MAP_vals[i], W2_sin_diff_plus, W2_sin_diff_minus, name=name))
					
		W1_sin_mxprbs = zeros(nwalkers)
		W2_sin_mxprbs = zeros(nwalkers)
					
		for i in range(nwalkers):
			W1_sin_mxprbs[i] = max(W1_sin_lnprobs[i])
			W2_sin_mxprbs[i] = max(W2_sin_lnprobs[i])

					
		print "W1 Max LnP = ", max(W1_sin_mxprbs)
		print "W2 Max LnP = ", max(W2_sin_mxprbs)		







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


		filename = "SHW1_results_%iwalkers.txt" %clen
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
Tt   =  (loadtxt("../dat/Lums_PG1302.dat", usecols=[0]))/(1.+zPG1302)
Lum  =  loadtxt("../dat/Lums_PG1302.dat", usecols=[1])
sigL =  loadtxt("../dat/Lums_PG1302.dat", usecols=[2])

Lum_Mean = mean(Lum)
Lum      = (Lum - Lum_Mean)  ## for def of mag - mag0

## put in time order
tLumS   = zip(Tt,Lum,sigL)
tLumS.sort()
TtLumS  = transpose(tLumS)
tsrt    =  TtLumS[0]
Lumsrt  =  TtLumS[1]
sigLsrt =  TtLumS[2]

tsrt = tsrt #- 49100
t_MJD = t_MJD #- 49100

Nt=40
ttopt = np.linspace(tsrt[0]-100, t_MJD[len(t_MJD)-1]+100,       Nt)


opti = -2.5*np.log10(Fsrc(ttopt*3600.*24, Dst, ma.pi/2., 0.0, Lav, betst, Inc, Ombn, 1.1)/FVbndRel)

ttopt = (ttopt*(1.+zPG1302) - 50000)
t_avg = (t_avg*(1.+zPG1302) - 50000)
tsrt  = (tsrt*(1.+zPG1302) - 50000)
t_MJD = (t_MJD*(1.+zPG1302) - 50000)

## z adjust time parameers
## PERIOD
#W1_sin_p_opt[1] = W1_sin_p_opt[1]*(1.+zPG1302) 
#W2_sin_p_opt[1] = W2_sin_p_opt[1]*(1.+zPG1302) 



plt.figure()
plt.errorbar(tsrt, Lumsrt+Lum_Mean-3.0, yerr=sigL, linestyle="none", color = "black", alpha=0.5) #alpha=0.1
Fs = plt.plot(ttopt, opti, linestyle = '--', color='blue', linewidth=2)


W1dat   = plt.errorbar(t_MJD, W1_mag, yerr=W1_sig, linestyle="none", color='orange', alpha=1., elinewidth=1.5)
W2dat   = plt.errorbar(t_MJD, W2_mag+0.5, yerr=W2_sig, linestyle="none", color='red', alpha=1., elinewidth=1.5)

W1av   = plt.errorbar(t_avg, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)
W2av   = plt.errorbar(t_avg, W2_avg+0.5, yerr=W2_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)


if (SinFit):
	W1sinsoln = plt.plot(ttopt, sinPoint(W1_sin_p_opt, (ttopt+50000)/(1.+zPG1302)), linestyle = '--', color='orange', linewidth=2)
	W2sinsoln = plt.plot(ttopt, sinPoint(W2_sin_p_opt, (ttopt+50000)/(1.+zPG1302))+0.5, linestyle = '--', color='red', linewidth=2)
if (ShellFit):
	if (fmin_Fit):
		W1shell = plt.plot(ttopt, magPoint_Shell(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_Shell(ShW2_p_opt, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
	else:
		W1shell = plt.plot(ttopt, magPoint_Shell(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_Shell(ShW2_p0_0, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
if (ThickFit):	
	if (fmin_Fit):
		W1shell = plt.plot(ttopt, magPoint_Thick(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_Thick(ShW2_p_opt, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
	else:
		W1shell = plt.plot(ttopt, magPoint_Thick(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
		W2shell = plt.plot(ttopt, magPoint_Thick(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
if (NoFit and pltShell):
	W1shell = plt.plot(ttopt, magPoint_Shell(ShW1_p0_0, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
	W2shell = plt.plot(ttopt, magPoint_Shell(ShW2_p0_0, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
if (NoFit and pltThick):
	W1shell = plt.plot(ttopt, magPoint_Thick(ShW1_p_opt, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
	W2shell = plt.plot(ttopt, magPoint_Thick(ShW2_p_opt, (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)
	


plt.grid(b=True, which='both')
		#plt.legend( [ s1[0], IR1[0], s2[0], IR2[0], s3[0], IR3[0]  ], (r'$i=0$','',   r'$i=\pi/4$','',   r'$i=\pi/2$', ''), loc='upper right')

plt.xlabel(r"$t$ [MJD]")
plt.ylabel("mag")
#plt.xlim(52000, 57500)
plt.xlim(3000, 8000)
#plt.ylim(10.5, 11.5)
#plt.ylim(plt.ylim(10.5, 12.3)[::-1])

		#plt.show()
plt.savefig("../emcee_data/"+Shell_File+"BestFit.png")
plt.clf()
