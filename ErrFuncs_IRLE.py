







##FOR MCMC
def magPoint_Shell(params, t, THEargs, RHStable, Ttable, rem_is_Rin):
	if (rem_is_Rin):
		sinJJ, cosTT, Rin, n0 = params
		rem_pls = 0.0
	else:
		sinJJ, cosTT, rem_pls, Rin, n0 = params

	n0 = n0 * 1.4032428247438431e-09
	Rin = Rin * 2.73213149e+18
	rem = Rin + rem_pls * 2.73213149e+18

	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
	
	FRel, numin, numax, Dist, Lav, Ombn, alph, pp, Rout,  aeff, nu0, nne, beta = THEargs
	IncFit = np.arccos(0.067/beta)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
	
	return -2.5*np.log10(Fobs_Shell(numin, numax, rem, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)


	##FOR MCMC
def magPoint_Thick_fmin(params, t, THEargs, RHStable, Ttable):
	#beta, cosJJ, Rin, thetT, n0 = params
	sinJJ, cosTT, Rin, pp, n0 = params
	if (sinJJ<-1.0 or sinJJ >1.0):
		return np.inf
	elif (cosTT<0.0 or cosTT >1.0):
		return np.inf
	else:
		Rin = Rin * 2.73213149e+18
		n0 = n0 * 1.4032428247438431e-09
		t = t * 86400.
		JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
		thetT = np.arccos(cosTT)
		
		#FRel, numin, numax, Dist, Lav, Ombn, alph, Rin, Rout, aeff, nu0, nne, beta = THEargs
		FRel, numin, numax, Dist, Lav, Ombn, alph, Rout, aeff, nu0, nne, beta = THEargs
		IncFit = np.arccos(0.067/beta)

		Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
		#return -2.5*np.log10(Fobs_Shell(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)
		return -2.5*np.log10(Fobs_Thick(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)


def magPoint_Thick(params, t, THEargs, RHStable, Ttable):
	#beta, cosJJ, Rin, thetT, n0 = params
	sinJJ, cosTT, Rin, pp, n0 = params
	Rin = Rin * 2.73213149e+18
	n0 = n0 * (pp-1.)/(ma.pi * aeff*aeff * Rde)   # must be greater than one
	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
		
	#FRel, numin, numax, Dist, Lav, Ombn, alph, Rin, Rout, aeff, nu0, nne, beta = THEargs
	FRel, numin, numax, Dist, Lav, Ombn, alph, Rout, aeff, nu0, nne, beta = THEargs
	IncFit = np.arccos(0.067/beta)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
	#return -2.5*np.log10(Fobs_Shell(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)
	return -2.5*np.log10(Fobs_Thick(numin, numax, t, Dist, Rout, Aargs, RHStable, Ttable)/FRel)
