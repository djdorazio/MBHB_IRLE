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

from scipy import optimize
from scipy.optimize import fmin

from FluxFuncs_IRLE import *
from ErrFuncs_IRLE import *
from emcee_Funcs import *

################################
###############################
### OPTIONS
################################
################################
Flat = False

SinFit = True
No_Prd = True

##multiprocessing
NThread = 4




zPG1302 = 0.2784






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



###### get average value for each cluster of data points in time
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
	Nseg1 = len(W1_sig[iseg[i]+1:iseg[i+1]])
	Nseg2 = len(W2_sig[iseg[i]+1:iseg[i+1]])

	#W1_avsg.append(np.sqrt(sum( (W1_sig[iseg[i]+1:iseg[i+1]])**2 ))/Nseg)
	#W2_avsg.append(np.sqrt(sum( (W2_sig[iseg[i]+1:iseg[i+1]])**2 ))/Nseg)

	W1_avsg.append(np.sqrt(sum( (W1_sig[iseg[i]+1:iseg[i+1]])**2 )/Nseg1  ))
	W2_avsg.append(np.sqrt(sum( (W2_sig[iseg[i]+1:iseg[i+1]])**2 )/Nseg2  ))

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





















if (Flat):
	ndim=1
	Shell_File = "FLAT_"
	param_names = ['mag0']
	Sinp0_W1 = [11.3]
	Sinp0_W2 = [10.3]
elif (SinFit):
	SinPrd = 1884.0/(1. + zPG1302)
	if (No_Prd):
		ndim = 3
		Shell_File = "SinFit_NOPrd"
		param_names = ['Amp','phase','mag0']
		Sinp0_src = [0.1269, 1403.4292, 14.8330]
		Sinp0_W1 = [0.0887, 0.01, 11.3]
		Sinp0_W2 = [0.1, 0.01, 10.3]	
	else:
		ndim = 4
		Shell_File = "SinFit_wPrd"
		param_names = ['Amp','Prd','phase','mag0']	
		Sinp0_src = [0.1269, SinPrd, 1403.4292, 14.8330]
		Sinp0_W1 = [0.0887, SinPrd, 0.01, 11.3]
		Sinp0_W2 = [0.1, SinPrd, 0.01, 10.3]
	
nwalkers = ndim*16

src_sin_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_Sinposterior, threads=NThread,args=(tsrt, Lumsrt, sigLsrt, Flat, SinFit, No_Prd))

W1_sin_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_Sinposterior, threads=NThread,args=(t_avg, W1_avg, W1_avsg, Flat, SinFit, No_Prd))
W2_sin_sampler  = emcee.EnsembleSampler(nwalkers, ndim, ln_Sinposterior, threads=NThread,args=(t_avg, W2_avg, W2_avsg, Flat, SinFit, No_Prd))


#p0 = np.array(p0)
src_sin_p0 = np.array(Sinp0_src)
W1_sin_p0 = np.array(Sinp0_W1)
W2_sin_p0 = np.array(Sinp0_W2)
src_sin_walker_p0 = np.random.normal(W1_sin_p0, np.abs(W1_sin_p0)*1E-4, size=(nwalkers, ndim))
W1_sin_walker_p0 = np.random.normal(W1_sin_p0, np.abs(W1_sin_p0)*1E-4, size=(nwalkers, ndim))
W2_sin_walker_p0 = np.random.normal(W2_sin_p0, np.abs(W2_sin_p0)*1E-4, size=(nwalkers, ndim))

			
clen = 4096
W1_sin_pos,_,_ = W1_sin_sampler.run_mcmc(W1_sin_walker_p0 , clen)

W2_sin_pos,_,_ = W2_sin_sampler.run_mcmc(W2_sin_walker_p0 , clen)

src_sin_pos,_,_ = src_sin_sampler.run_mcmc(src_sin_walker_p0 , clen)

print "SAVING THE PICKLE mmmmm"
with open("../emcee_data/Pickles/PG1302_src_sin_%iwalkers.pickle" %clen, "w") as f1:
	pickle.dump((src_sin_sampler.chain, src_sin_sampler.lnprobability), f1)

with open("../emcee_data/Pickles/PG1302_W1_sin_%iwalkers.pickle" %clen, "w") as f2:
	pickle.dump((W1_sin_sampler.chain, W1_sin_sampler.lnprobability), f2)

with open("../emcee_data/Pickles/PG1302_W2_sin_%iwalkers.pickle" %clen, "w") as f3:
	pickle.dump((W2_sin_sampler.chain, W2_sin_sampler.lnprobability), f3)



		


### OPEN OUTPUT DATA
with open("../emcee_data/Pickles/PG1302_src_sin_%iwalkers.pickle" %clen) as f1:
	src_sin_chain,src_sin_lnprobs = pickle.load(f1)

with open("../emcee_data/Pickles/PG1302_W1_sin_%iwalkers.pickle" %clen) as f2:
	W1_sin_chain,W1_sin_lnprobs = pickle.load(f2)

with open("../emcee_data/Pickles/PG1302_W2_sin_%iwalkers.pickle" %clen) as f3:
	W2_sin_chain,W2_sin_lnprobs = pickle.load(f3)




src_sin_flatchain = np.vstack(src_sin_chain[:,clen/4:])
src_sin_flatlnprobs = np.vstack(src_sin_lnprobs[:,clen/4:])

W1_sin_flatchain = np.vstack(W1_sin_chain[:,clen/4:])
W1_sin_flatlnprobs = np.vstack(W1_sin_lnprobs[:,clen/4:])

W2_sin_flatchain = np.vstack(W2_sin_chain[:,clen/4:])
W2_sin_flatlnprobs = np.vstack(W2_sin_lnprobs[:,clen/4:])

			

src_sin_p_opt  = src_sin_flatchain[src_sin_flatlnprobs.argmax()]		
W1_sin_p_opt  = W1_sin_flatchain[W1_sin_flatlnprobs.argmax()]
W2_sin_p_opt  = W2_sin_flatchain[W2_sin_flatlnprobs.argmax()]













print "ANALYSING MCMC (!)..."
with open("../emcee_data/Pickles/PG1302_src_sin_%iwalkers.pickle" %clen) as f1:
	src_sin_chain,src_sin_lnprobs = pickle.load(f1)
with open("../emcee_data/Pickles/PG1302_W1_sin_%iwalkers.pickle" %clen) as f2:
	W1_sin_chain,W1_sin_lnprobs = pickle.load(f2)
with open("../emcee_data/Pickles/PG1302_W2_sin_%iwalkers.pickle" %clen) as f3:
	W2_sin_chain,W2_sin_lnprobs = pickle.load(f3)


##PLOT dem WALKERS
for k in range(src_sin_chain.shape[2]):
	plt.figure()
	#plt.figure()
	for i in range(src_sin_chain.shape[0]):
		plt.plot(src_sin_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
		plt.ylabel(param_names[k])
		plt.xlabel('steps')
	plt.savefig('../emcee_data/src_'+Shell_File+'_PG1302_%s_%iwalkers.png' %(param_names[k],clen))
	plt.clf()

for k in range(W1_sin_chain.shape[2]):
	plt.figure()
	#plt.figure()
	for i in range(W1_sin_chain.shape[0]):
		plt.plot(W1_sin_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
		plt.ylabel(param_names[k])
		plt.xlabel('steps')
	plt.savefig('../emcee_data/W1_'+Shell_File+'_PG1302_%s_%iwalkers.png' %(param_names[k],clen))
	plt.clf()

for k in range(W2_sin_chain.shape[2]):
	plt.figure()
	#plt.figure()
	for i in range(W2_sin_chain.shape[0]):
		plt.plot(W2_sin_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
		plt.ylabel(param_names[k])
		plt.xlabel('steps')
	plt.savefig('../emcee_data/W2_'+Shell_File+'_PG1302_%s_%iwalkers.png' %(param_names[k],clen))
	plt.clf()



###CORNER PLOT	
src_sin_flatchain = np.vstack(src_sin_chain[:,clen/2:])
src_sin_flatlnprobs = np.vstack(src_sin_lnprobs[:,clen/2:])

W1_sin_flatchain = np.vstack(W1_sin_chain[:,clen/2:])
W1_sin_flatlnprobs = np.vstack(W1_sin_lnprobs[:,clen/2:])

W2_sin_flatchain = np.vstack(W2_sin_chain[:,clen/2:])
W2_sin_flatlnprobs = np.vstack(W2_sin_lnprobs[:,clen/2:])

		

#import triangle
import corner as triangle
src_sin_fig = triangle.corner(src_sin_flatchain, labels=param_names)			
src_sin_fig.savefig('../emcee_data/src_'+Shell_File+'_PG1302_Corner_Plot_%iwalkers.png' %clen)


W1_sin_fig = triangle.corner(W1_sin_flatchain, labels=param_names)			
W1_sin_fig.savefig('../emcee_data/W1_'+Shell_File+'_PG1302_Corner_Plot_%iwalkers.png' %clen)


W2_sin_fig = triangle.corner(W2_sin_flatchain, labels=param_names)			
W2_sin_fig.savefig('../emcee_data/W2_'+Shell_File+'_PG1302Corner_Plot_%iwalkers.png' %clen)



## Do some stats on the walkers
from scipy.stats import scoreatpercentile as scoretpercentile



## max posterior + percentiles
src_sin_MAP_vals = src_sin_flatchain[src_sin_flatlnprobs.argmax()]
src_sin_perc = scoretpercentile(src_sin_flatchain, [15,85], axis=0)

W1_sin_MAP_vals = W1_sin_flatchain[W1_sin_flatlnprobs.argmax()]
W1_sin_perc = scoretpercentile(W1_sin_flatchain, [15,85], axis=0)

W2_sin_MAP_vals = W2_sin_flatchain[W2_sin_flatlnprobs.argmax()]
W2_sin_perc = scoretpercentile(W2_sin_flatchain, [15,85], axis=0)



filename = "Sin_results_"+Shell_File+"%iwalkers.txt" %clen
print "Printing Results"
target = open(filename, 'w')
target.truncate()


for i,name in enumerate(param_names):
	src_sin_diff_minus = src_sin_MAP_vals[i] - src_sin_perc[0,i]
	src_sin_diff_plus = src_sin_perc[1,i] - src_sin_MAP_vals[i]
	target.write("Src: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(src_sin_MAP_vals[i], src_sin_diff_plus, src_sin_diff_minus, name=name))
	target.write("\n")

target.write("\n")
for i,name in enumerate(param_names):
	W1_sin_diff_minus = W1_sin_MAP_vals[i] - W1_sin_perc[0,i]
	W1_sin_diff_plus = W1_sin_perc[1,i] - W1_sin_MAP_vals[i]
	target.write("W1: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(W1_sin_MAP_vals[i], W1_sin_diff_plus, W1_sin_diff_minus, name=name))
	target.write("\n")

target.write("\n")
for i,name in enumerate(param_names):
	W2_sin_diff_minus = W2_sin_MAP_vals[i] - W2_sin_perc[0,i]
	W2_sin_diff_plus = W2_sin_perc[1,i] - W2_sin_MAP_vals[i]
	target.write("W2: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(W2_sin_MAP_vals[i], W2_sin_diff_plus, W2_sin_diff_minus, name=name))
	target.write("\n")	


if (No_Prd):
	target.write("Period = %4g" %SinPrd)	
	target.write("\n")		


src_sin_mxprbs = zeros(nwalkers)
W1_sin_mxprbs = zeros(nwalkers)
W2_sin_mxprbs = zeros(nwalkers)
			

for i in range(nwalkers):
	src_sin_mxprbs[i] = max(src_sin_lnprobs[i])
	W1_sin_mxprbs[i] = max(W1_sin_lnprobs[i])
	W2_sin_mxprbs[i] = max(W2_sin_lnprobs[i])

chi2_pdf_src = -max(src_sin_mxprbs)/(len(Lumsrt) - len(param_names) - 1)
chi2_pdf_W1 = -max(W1_sin_mxprbs)/(len(W1_avg) - len(param_names) - 1)
chi2_pdf_W2 = -max(W2_sin_mxprbs)/(len(W2_avg) - len(param_names) - 1)
target.write("\n")		
target.write("Shell src reduced chi2 =  %04g" %chi2_pdf_src)
target.write("\n")		
target.write("Shell W1 reduced chi2 =  %04g" %chi2_pdf_W1)
target.write("\n")
target.write("Shell W2 reduced chi2 =  %04g" %chi2_pdf_W2)
target.write("\n")
target.write("\n")


AV = src_sin_p_opt[0]
AW1 = W1_sin_p_opt[0]
AW2 = W2_sin_p_opt[0]


AW1oAV = AW1/AV
AW2oAV = AW2/AV

dAV_dw = (src_sin_MAP_vals[0] - src_sin_perc[0,0])
dAV_up = (src_sin_perc[1,0] - src_sin_MAP_vals[0])

dW1_dw = (W1_sin_MAP_vals[0] - W1_sin_perc[0,0])
dW1_up = (W1_sin_perc[1,0] - W1_sin_MAP_vals[0])

dW2_dw = (W2_sin_MAP_vals[0] - W2_sin_perc[0,0])
dW2_up = (W2_sin_perc[1,0] - W2_sin_MAP_vals[0]) 

delAV = 0.5* ( (src_sin_MAP_vals[0] - src_sin_perc[0,0]) + (src_sin_perc[1,0] - src_sin_MAP_vals[0]) )
			
delAW1 = 0.5* ( (W1_sin_MAP_vals[0] - W1_sin_perc[0,0]) + (W1_sin_perc[1,0] - W1_sin_MAP_vals[0]) )
	
delAW2 = 0.5* ( (W2_sin_MAP_vals[0] - W2_sin_perc[0,0]) + (W2_sin_perc[1,0] - W2_sin_MAP_vals[0]) )


delW1oV = ma.sqrt( (delAW1/AV)**2 + (AW1/AV/AV * delAV)**2 )
delW1oV_up = ma.sqrt( (dW1_up/AV)**2 + (AW1/AV/AV * dAV_up )**2 )
delW1oV_dw = ma.sqrt( (dW1_dw/AV)**2 + (AW1/AV/AV * dAV_dw )**2 )

delW2oV = ma.sqrt( (delAW2/AV)**2 + (AW2/AV/AV * delAV)**2 )
delW2oV_up = ma.sqrt( (dW2_up/AV)**2 + (AW2/AV/AV * dAV_up)**2 )
delW2oV_dw = ma.sqrt( (dW2_dw/AV)**2 + (AW2/AV/AV * dAV_dw)**2 )

target.write("AW1/AV = %g p %g  m %g" %(AW1oAV, delW1oV_up, delW1oV_dw))
target.write("\n")
target.write("AW2/AV = %g p %g  m %g" %(AW2oAV, delW2oV_up, delW2oV_dw))

print "AW1/AV = %g p %g  m %g" %(AW1oAV, delW1oV_up, delW1oV_dw)
print "AW2/AV = %g p %g  m %g" %(AW2oAV, delW2oV_up, delW2oV_dw)

target.close()


# print "PLOTTING BEST FIT LIGHT CURVES"
from Gen_Plot import *
# if (Shell_OptThin):
# 	#Plot_Shell_Thin_ISO(p_opt, 40, Shell_File,   W1args, W2args, RHS_table, T_table,  tsrt, t_avg, t_MJD,    Lumsrt, W1_mag, W2_mag, W1_avg, W2_avg,   sigL, W1_sig, W2_sig, W1_avsg, W2_avsg)
# 	Plot_SinFlat(W1_sin_p_opt, W1_sin_p_opt, 40, Shell_File,   tsrt, t_avg, t_MJD,    Lumsrt, W1_mag, W2_mag, W1_avg, W2_avg,   sigL, W1_sig, W2_sig, W1_avsg, W2_avsg)

print "Plotting Best Fit Sins"
#p_Fsrc_sin = [0.1269, 1403.4292, 14.8330]
Nt = 40

#Plot_Sin(p_Fsrc_sin, W1_sin_p_opt, W2_sin_p_opt, Nt, Shell_File,    Flat, SinFit, No_Prd,     tsrt, t_avg, t_MJD,    Lumsrt, W1_mag, W2_mag, W1_avg, W2_avg,   sigL, W1_sig, W2_sig ,W1_avsg, W2_avsg)
Plot_Sin(src_sin_p_opt, W1_sin_p_opt, W2_sin_p_opt, Nt, Shell_File,    Flat, SinFit, No_Prd,     tsrt, t_avg, t_MJD,    Lumsrt, W1_mag, W2_mag, W1_avg, W2_avg,   sigL, W1_sig, W2_sig ,W1_avsg, W2_avsg)















