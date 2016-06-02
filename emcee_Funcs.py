import numpy as np
from numpy import *
from ErrFuncs_IRLE import *

### PRIORS
def ln_prior(params):
	sinJJ, cosTT, Rin, alpha = params
					
	if sinJJ < -1 or sinJJ > 1:
		return -np.inf

	if cosTT < 0 or cosTT > 1:
		return -np.inf
						
	if Rin <= 0.0:
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

##LIKLIHOODS
def ln_SHThin_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			return -(OpThin_TorShell_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		
##LIKLIHOODS
def ln_ISO_SHThin_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			return -(ISO_OpThin_TorShell_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2))
		


def ln_Sinlikelihood(p, t, y, dy):
			return -(SinErr2(p, t, y, dy)) 







##POSTERIORS

def ln_SHThin_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			ln_p = ln_prior(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_SHThin_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
			return ln_l + ln_p


def ln_ISO_SHThin_posterior(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
			ln_p = ln_prior(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_ISO_SHThin_likelihood(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2)
			return ln_l + ln_p



def ln_Sinposterior(p, t, y, dy):
			ln_p = ln_Sinprior(p)
			if not np.isfinite(ln_p):
				return -np.inf
			
			ln_l = ln_Sinlikelihood(p, t, y, dy)
			return ln_l + ln_p