import cPickle as pickle


import matplotlib
matplotlib.use('Agg')

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt


matplotlib.rcParams.update({'font.size': 16})


from FluxFuncs_IRLE import *
from ErrFuncs_IRLE import *


import scipy.integrate as intgt

from scipy import special as spc

NThread = 8
AddW4 = False #W1, W2, W3, W4 if True, otherwise see AddW3
AddW3 = True #W1, W2, W3 if True, #W1, W2 if False


##GLOBAL PHYSICS CONSTANTS (cgs):
c = 2.9979*10**(10)
sigSB = 5.670*10**(-5)

##ISO VARS:
Amp = 0.15# For alph = 1, bet = 0.07 0.22105  ##For alpha=0.0, beta=0.0677
t0 = ma.pi/2. #0.0




#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
Omb = 1./(1*yr2sec)
L0 = 6.78*10**46
MPGmx = 10**9.4*Msun
Ryr = c*yr2sec
RdPG = 1.0*pc2cm#ma.sqrt(1.)*2.8 *pc2cm
OmPG = 2.*ma.pi/(1474.*3600*24) #Omb*2.*ma.pi/4.1
Ombn = OmPG
#alphnu = 1.1

Rorb = c*2.*ma.pi/Omb
Ompc = 2.*ma.pi*c/pc2cm/2.


## TEST VALUES
### DUST stuff
## for Qv
nne = 0.0#1.8
nu0 = numicron#/0.37

Rde = RdPG
Rrout = 1.0*Rde
pp = 2.0
thetTst = 0.0#1.*ma.pi/2.


aeff = (c/nu0)/(2.*ma.pi)
#aeff = 0.16*10**(-4) #(0.1 micrometer is an average ISM dust grain size - choose 0.16 to make nu0~1um)
#md = 10**(-14)
n10 = 1.0/(ma.pi*Rde*aeff*aeff) * (pp-1.)
nfac = 10.0
n0 = nfac*n10 ##6.*10**5*Msun/md * 1./(4./3.*ma.pi*(Rrout**3 - Rde**3))


Lav = L0
betst = 0.07  ## gives 0.14 mag amplitdue at I=0, alpha=1.1
Inc = 0.0#ma.pi/2#ma.acos(0.07/betst)#0.*np.pi/4.
#Ombn = 2.*ma.pi/(Rde/c)  ##(2pi/P)
alph = 1.0  ## for getting optical LC
Dopalph = -1.0 ## for illuminate dust
Dst = 1.4*10**9*pc2cm




### WISE BAND + Observational STUFF
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3


W1_mid = (W1mx + W1mn)/2.
W2_mid = (W2mx + W2mn)/2.

W3_mid = numicron/12.
W4_mid = numicron/22.



nuVbnd = c/(5.45*10**(-5))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-21)*8.8560*10**(13)#(W1mn + W1mx)/2
FW2Rel = 1.71787*10**(-21)*6.4451*10**(13)#(W2mn + W2mx)/2


# if (not (Qv_k or Qv_nu0)):
# 	##TABULATE T's and RHSs
# 	print "Creating look up tables"
# 	NT = 10000
# 	RHS_table = np.zeros(NT)
# 	T_table = np.linspace(1., 2000., NT)
# 	for i in range(NT):
# 		RHS_table[i] = T_RHS(T_table[i], nu0, nne)







Tmin = 100.
Tsub = 1800.
NT   = 1800

















####MAIN




#read in fluxes and errors, take average
FnuW1_t = np.genfromtxt("../dat/W1.txt",usecols=0) * 10.**(-26)
FnuW2_t = np.genfromtxt("../dat/W2.txt",usecols=0) * 10.**(-26)
ErrW1_t = np.genfromtxt("../dat/W1.txt",usecols=1) * 10.**(-26)
ErrW2_t = np.genfromtxt("../dat/W2.txt",usecols=1) * 10.**(-26)


Fw1 = np.mean(FnuW1_t) 
Fw2 = np.mean(FnuW2_t)
Fw3 = 36.0745 * 10.**(-26)
Fw4 = 101.8591 * 10.**(-26)


ErrW1 = np.mean(ErrW1_t)
ErrW2 = np.mean(ErrW2_t)
ErrW3 = 0.5647 * 10.**(-26)
ErrW4 = 2.6262 * 10.**(-26)


## BEST FIT Black Body
from scipy import optimize
from scipy.optimize import fmin

from scipy import special as spc
from emcee_Funcs import *

if (AddW4):
	# All
	Shell_File = "AddW4"
	nus = (W1_mid, W2_mid, W3_mid, W4_mid)
	nus = np.array(nus)
	Flxs = (Fw1, Fw2, Fw3, Fw4)
	Flxs = np.array(Flxs)
	Errs = (ErrW1, ErrW2, ErrW3, ErrW4)
	Errs = np.array(Errs)
elif (AddW3):
	# NO W4
	Shell_File = "AddW3"
	nus = (W1_mid, W2_mid, W3_mid)
	nus = np.array(nus)
	Flxs = (Fw1, Fw2, Fw3)
	Flxs = np.array(Flxs)
	Errs = (ErrW1, ErrW2, ErrW3)
	Errs = np.array(Errs)
else:
	#No W# and no W4
	Shell_File = "W1W2"
	nus = (W1_mid, W2_mid)
	nus = np.array(nus)
	Flxs = (Fw1, Fw2)
	Flxs = np.array(Flxs)
	Errs = (ErrW1, ErrW2)
	Errs = np.array(Errs)	





import emcee

if (AddW4):
	Shell_File = "W1W2W3W4_2_param_BBfit_"
if (AddW3):
	Shell_File = "W1W2W3_2_param_BBfit_"
else:
	Shell_File = "W1W2_2_param_BBfit_"





param_names = [r"$T_d$", r"$X$"]

XX0 = 4.*ma.pi ## XX=4 pi fcov (R/D)^2, mult by (1.e-9)**2 in ERrFunc
BB_p0 = [577.0, XX0]

ndim = len(BB_p0)
nwalkers = ndim*16#32


BB_sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_BBCorposterior, threads=NThread, args=(nus, Flxs, Errs) )



BB_p0 = np.array(BB_p0)
BB_walker_p0 = np.random.normal(BB_p0, np.abs(BB_p0)*1E-4, size=(nwalkers, ndim))

clen = 2#4096
BB_pos,_,_ = BB_sampler.run_mcmc(BB_walker_p0 , clen)


print "SAVING THE PICKLE mmmmm"
with open("../emcee_data/Pickles/PG1302_BBCorfit_%iwalkers.pickle"%clen+Shell_File, "w") as f1:
	pickle.dump((BB_sampler.chain, BB_sampler.lnprobability), f1)

### OPEN OUTPUT DATA
with open("../emcee_data/Pickles/PG1302_BBCorfit_%iwalkers.pickle"%clen+Shell_File) as f1:
	BB_chain,BB_lnprobs = pickle.load(f1)


BB_flatchain = np.vstack(BB_chain[:,clen/4:])
BB_flatlnprobs = np.vstack(BB_lnprobs[:,clen/4:])

			
BB_p_opt  = BB_flatchain[BB_flatlnprobs.argmax()]

print "ANALYSING MCMC (!)..."
with open("../emcee_data/Pickles/PG1302_BBCorfit_%iwalkers.pickle"%clen+Shell_File) as f1:
	BB_chain,BB_lnprobs = pickle.load(f1)


##PLOT dem WALKERS
for k in range(BB_chain.shape[2]):
	plt.figure()
	for i in range(BB_chain.shape[0]):
		plt.plot(BB_chain[i,:,k], drawstyle='steps', color='k', marker=None, alpha=0.2)
		plt.ylabel(param_names[k])
		plt.xlabel('steps')
		#print param_names[k]
		plt.tight_layout()
	plt.savefig('../emcee_data/BBCorfit_PG1302_%s_%iwalkers'%(param_names[k],clen)+Shell_File+'.png')
	plt.clf()


###CORNER PLOT	
BB_flatchain = np.vstack(BB_chain[:,clen/4:])
BB_flatlnprobs = np.vstack(BB_lnprobs[:,clen/4:])

#import triangle
import corner as triangle
BB_fig = triangle.corner(BB_flatchain, labels=param_names, quantiles=[0.15, 0.5, 0.85],show_titles=True, title_kwargs={"fontsize": 14},label_kwargs={"fontsize": 20})	
BB_fig.savefig('../emcee_data/BBCorfit_PG1302_Corner_Plot_%iwalkers'%clen+Shell_File+'.png')


## Do some stats on the walkers
from scipy.stats import scoreatpercentile as scoretpercentile
## max posterior + percentiles
BB_MAP_vals = BB_flatchain[BB_flatlnprobs.argmax()]
BB_perc = scoretpercentile(BB_flatchain, [15,85], axis=0)


filename = "BBfit_resuts_%iwalkers"%clen+Shell_File+".txt" 
print "Printing Results"
target = open(filename, 'w')
target.truncate()


for i,name in enumerate(param_names):
	BB_diff_minus = BB_MAP_vals[i] - BB_perc[0,i]
	BB_diff_plus = BB_perc[1,i] - BB_MAP_vals[i]
	target.write("BBfit: {name}: {0:.4f} + {1:.4f} - {2:.4f}".format(BB_MAP_vals[i], BB_diff_plus, BB_diff_minus, name=name))
	target.write("\n")





BB_mxprbs = zeros(nwalkers)

			

for i in range(nwalkers):
	BB_mxprbs[i] = max(BB_lnprobs[i])


chi2_pdf_BB = -max(BB_mxprbs)#/(len(nus) - len(param_names) - 1)

target.write("\n")		
target.write("BBfit reduced chi2 =  %04g" %chi2_pdf_BB)

target.close()




### GET Rd, cos(thT), BC from XX and Td fit

### get avg deltas
delT   = 0.5* ( (BB_MAP_vals[0] - BB_perc[0,0]) + (BB_perc[1,0] - BB_MAP_vals[0]) )
delXX = 0.5* ( (BB_MAP_vals[1] - BB_perc[0,1]) + (BB_perc[1,1] - BB_MAP_vals[1]) ) *(1.e-9)**2*(1000.)**4



Td = BB_p_opt[0]
XX = (BB_p_opt[1]+delXX) *(1.e-9)**2*(1000.)**4




BCmean = 0.9
BCvar = 0.3
BCmin = BCmean - BCvar
BCmax = BCmean + BCvar


fmean = 4.*XX * sigSB /(BCmean*Lav) * Dst*Dst 
fmin = 4.*XX * sigSB /(BCmax*Lav) * Dst*Dst 
fmax = 4.*XX * sigSB /(BCmin*Lav) * Dst*Dst 
delfmx = fmax - fmean
delfmn = fmean - fmin

Rdmean = np.sqrt(XX/Td**4/(4.*ma.pi * fmean))* Dst/pc2cm
Rdmin = np.sqrt(XX/Td**4/(4.*ma.pi * fmax)) * Dst/pc2cm
Rdmax = np.sqrt(XX/Td**4/(4.*ma.pi * fmin)) * Dst/pc2cm
delRmx = Rdmax - Rdmean
delRmn = Rdmean - Rdmin

tdop_mean = Rdmean*pc2cm/c/(2.*ma.pi) * OmPG
tdop_min = Rdmin*pc2cm/c/(2.*ma.pi) * OmPG
tdop_max = Rdmax*pc2cm/c/(2.*ma.pi) * OmPG
delTdoPmx = tdop_max - tdop_mean
delTdoPmn = tdop_mean - tdop_min



### WRITE RESULTS TO FILE
filename = "../emcee_Results/BBCorfit_BC%gpm%g_%iwalkers"%(BCmean, BCvar, clen)+Shell_File+"_FULL_2Dparam_RESULTS.txt" 
target = open(filename, 'w')
target.truncate()
target.write(Shell_File)
target.write("\n")	
target.write("T_opt = %g +- %g K" %(Td, delT))
target.write("\n")	
target.write("XX = %g +- %g" % (XX, delXX))
target.write("\n")	
target.write("Fcover = %g - %g + %g" %(fmean, delfmn, delfmx))
target.write("\n")	
target.write("R_opt = %g - %g + %g pc"  %(Rdmean, delRmn, delRmx) )
target.write("\n")	
target.write( "t_d/P = %g - %g + %g"  %(tdop_mean, delTdoPmn, delTdoPmx) )
target.close()





########PLOTTING

nu = np.linspace(0.01*numicron, 1*numicron, 100)/10**14



plt.figure()
plt.scatter(W1_mid/10**14, Fw1* 10.**(26), color='orange', s=40, marker='o')
plt.scatter(W2_mid/10**14, Fw2* 10.**(26), color='red', s =40, marker='o')
plt.errorbar(W1_mid/10**14, Fw1* 10.**(26), yerr=ErrW1*10.**(26), linestyle="none", color='orange', alpha=1., elinewidth=1.5)
plt.errorbar(W2_mid/10**14, Fw2* 10.**(26), yerr=ErrW2*10.**(26), linestyle="none", color='red', alpha=1., elinewidth=1.5)

if (AddW4):
	plt.scatter(W3_mid/10**14, Fw3* 10.**(26), color='purple', s=40, marker='o')
	plt.scatter(W4_mid/10**14, Fw4* 10.**(26), color='brown', s =40, marker='o', alpha=1.0)
	plt.errorbar(W3_mid/10**14, Fw3* 10.**(26), yerr=ErrW3*10.**(26), linestyle="none", color='purple', alpha=1., elinewidth=1.5)
	plt.errorbar(W4_mid/10**14, Fw4* 10.**(26), yerr=ErrW4*10.**(26), linestyle="none", color='brown', alpha=1.0, elinewidth=1.5)
elif(AddW3):
	plt.scatter(W3_mid/10**14, Fw3* 10.**(26), color='purple', s=40, marker='o')
	plt.scatter(W4_mid/10**14, Fw4* 10.**(26), color='brown', s =40, marker='o', alpha=0.5)
	plt.errorbar(W3_mid/10**14, Fw3* 10.**(26), yerr=ErrW3*10.**(26), linestyle="none", color='purple', alpha=1., elinewidth=1.5)
	plt.errorbar(W4_mid/10**14, Fw4* 10.**(26), yerr=ErrW4*10.**(26), linestyle="none", color='brown', alpha=0.5, elinewidth=1.5)
else:
	plt.scatter(W3_mid/10**14, Fw3* 10.**(26), color='purple', s=40, marker='o', alpha=0.5)
	plt.scatter(W4_mid/10**14, Fw4* 10.**(26), color='brown', s =40, marker='o', alpha=0.5)
	plt.errorbar(W3_mid/10**14, Fw3* 10.**(26), yerr=ErrW3*10.**(26), linestyle="none", color='purple', alpha=0.5, elinewidth=1.5)
	plt.errorbar(W4_mid/10**14, Fw4* 10.**(26), yerr=ErrW4*10.**(26), linestyle="none", color='brown', alpha=0.5, elinewidth=1.5)



plt.plot(nu,  Bv(nu*10**14, Td) * XX/Td**4 * 10.**(26), color = 'gray', linewidth = 2)

plt.xlabel(r'$\nu$ [$10^{14}$ Hz]')
plt.ylabel('Flux [mJy]')

plt.xlim(0.0,1.5)
plt.ylim(0.0,110.0)

plt.tight_layout()


if (AddW4):
	plt.savefig("../emcee_data/BBCorfit_BestFit_ModBlackBody_clen%g_%gwalkers_AddW4.png" %(clen, nwalkers))
if (AddW3):
	plt.savefig("../emcee_data/BBCorfit_BestFit_ModBlackBody_clen%g_%gwalkers_AddW3.png" %(clen, nwalkers))
else:
	plt.savefig("../emcee_data/BBCorfit_BestFit_ModBlackBody_clen%g_%gwalkers_W1W2.png" %(clen, nwalkers))


















