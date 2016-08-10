import numpy as np
from numpy import *
from ErrFuncs_IRLE import *

### PRIORS
# def ln_prior(params):
# 	#sinJJ, cosTT, Rin, alpha = params
# 	sinJJ, cosTT, Rin = params
					
# 	if sinJJ < -1 or sinJJ > 1:
# 		return -np.inf

# 	if cosTT < 0 or cosTT > 1:
# 		return -np.inf
						
# 	if Rin <= 0.0:
# 		return -np.inf
			
# 	return 0.

#PRIORS FROM MEASUREMENTS!
def ln_prior(params):
	#sinJJ, cosTT, Rin, alpha = params
	sinJJ, cosTT, Rin = params
					
	if sinJJ < -1.0 or sinJJ > 1.0:
		return -np.inf

	if cosTT < -1.0 or cosTT > 1.0:
		return -np.inf
						
	if Rin > 5.2 or Rin < 1.0 :
		return -np.inf
			
	return 0.


def ln_prior_mag0(params):
	sinJJ, cosTT, Rin, mag0_W1 = params
					
	if sinJJ < -1.0 or sinJJ > 1.0:
		return -np.inf

	if cosTT < -1.0 or cosTT > 1.0:
		return -np.inf
						
	if Rin > 5.2 or Rin < 1.0 :
		return -np.inf
			
	return 0.


def ln_prior_TwoRs(params):
	#sinJJ, cosTT, Rin, alpha = params
	sinJJ, cosTT, RW1, RW2 = params
					
	if sinJJ < -1.0 or sinJJ > 1.0:
		return -np.inf

	if cosTT < -1.0 or cosTT > 1.0:
		return -np.inf
						
	# if RW1 > 5.2 or RW1 < 1.0 :
	# 	return -np.inf

	# if RW2 > 5.2 or RW2 < 1.0 :
	# 	return -np.inf

	if RW1 < 0.0:
		return -np.inf

	if RW2 < 0.0 :
		return -np.inf
			
	return 0.


def ln_Sinprior(p, Flat, SinFit, No_Prd):
	if (Flat):
		return 0
	elif (SinFit):
		if (No_Prd):
			Prd = 1884./(1+0.2784)
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



def ln_BBprior(params):
	#sinJJ, cosTT, Rin, alpha = params
	Td, sqtfR  = params
					
	if Td < 0.0:
		return -np.inf

	if sqtfR < 0.0:
		return -np.inf
			
	return 0.



def ln_BBprior_Qv(params):
	#sinJJ, cosTT, Rin, alpha = params
	Td, nu0, gam, sqtfR  = params
	#Td, nu0, gam, fcov  = params
					
	if Td < 0.0:
		return -np.inf

	if nu0<=0.0 or nu0>10:
		return -np.inf
						
	if gam < 0.0 or gam>10.:
		return -np.inf

	if sqtfR < 0.0:
		return -np.inf
			
	#if fcov < 0.0 or fcov > 1.:
	#	return -np.inf


	return 0.


##LIKLIHOODS
def ln_SHThin_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			return -(OpThin_TorShell_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		

def ln_SHThin_likelihood_mag0(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			return -(OpThin_TorShell_Err2_mag0(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		

def ln_SHThin_likelihood_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			return -(OpThin_TorShell_Err2_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		

##LIKLIHOODS
def ln_ISO_SHThin_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			return -(ISO_OpThin_TorShell_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		
def ln_ISO_SHThin_likelihood_mag0(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			return -(ISO_OpThin_TorShell_Err2_mag0(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		

def ln_ISO_SHThin_likelihood_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			return -(ISO_OpThin_TorShell_Err2_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		

def sinPoint(params, t,  Flat, SinFit, No_Prd):
	if (Flat):
		return params
	elif (SinFit):
		if (No_Prd):
			Prd = 1884./(1.+0.2784)
			Amp, phs, mag0 = params
		else:
			Amp, Prd, phs, mag0 = params
		#Amp, phs, mag0 = params
		#Prd=1884.
	return Amp*np.sin( 2.*ma.pi/Prd*(t - phs) - np.pi/2) + mag0 

def SinErr2(p, t, y, dy, Flat, SinFit, No_Prd):
	print "EVAL", p
	chi = (y - sinPoint(p, t, Flat, SinFit, No_Prd) )/ dy
	chi2 = sum(chi*chi)
	print(chi2)
	return chi2

def ln_Sinlikelihood(p, t, y, dy, Flat, SinFit, No_Prd):
	return -(SinErr2(p, t, y, dy, Flat, SinFit, No_Prd)) 



def ln_BBlikelihood(p, t, y, dy):
	return -(BB_Err2(p, t, y, dy)) 


def ln_BBlikelihood_Qv(p, t, y, dy):
	return -(BB_Err2_Qv(p, t, y, dy)) 






##POSTERIORS

def ln_SHThin_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			ln_p = ln_prior(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_SHThin_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
			return ln_l + ln_p


def ln_SHThin_posterior_mag0(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			ln_p = ln_prior_mag0(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_SHThin_likelihood_mag0(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
			return ln_l + ln_p


def ln_SHThin_posterior_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			ln_p = ln_prior_TwoRs(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_SHThin_likelihood_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
			return ln_l + ln_p

def ln_ISO_SHThin_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			ln_p = ln_prior(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_ISO_SHThin_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
			return ln_l + ln_p

def ln_ISO_SHThin_posterior_mag0(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			ln_p = ln_prior_mag0(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_ISO_SHThin_likelihood_mag0(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
			return ln_l + ln_p


def ln_ISO_SHThin_posterior_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			ln_p = ln_prior_TwoRs(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_ISO_SHThin_likelihood_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
			return ln_l + ln_p



def ln_Sinposterior(p, t, y, dy, Flat, SinFit, No_Prd):
			ln_p = ln_Sinprior(p, Flat, SinFit, No_Prd)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_Sinlikelihood(p, t, y, dy, Flat, SinFit, No_Prd)
			return ln_l + ln_p


def ln_BBposterior(p, t, y, dy):
			ln_p = ln_BBprior(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_BBlikelihood(p, t, y, dy)
			return ln_l + ln_p


def ln_BBposterior_Qv(p, t, y, dy):
			ln_p = ln_BBprior_Qv(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_BBlikelihood_Qv(p, t, y, dy)
			return ln_l + ln_p





