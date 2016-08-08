import numpy as np
import math as ma
#import numexpr as ne
import scipy as sc
from FluxFuncs_IRLE import *
import time




def BB_Err2(p, nu, y, dy):
	Td, sqtfR  = p

	gam = 1.8

	nu0 = numicron/0.37

	sqtfR = sqtfR#*pc2cm
	Dst = 1.4*10**9#*pc2cm
	Lav = 6.78*10**46
	print p

	pref = np.ones(len(nu))
	for i in range(len(nu)):
		pref[i] = min(1., (nu[i]/nu0)**(gam))
	chi = (y - pref*Bv(nu, Td)* (sqtfR/Dst)**2 )/ dy
	chi2 = sum(chi*chi)
	print chi2
	return chi2



def BB_Err2_Qv(p, nu, y, dy):
	Td, nu0, gam, sqtfR  = p

	sqtfR = sqtfR#*pc2cm
	nu0 = nu0*10**14
	Dst = 1.4*10**9#*pc2cm
	Lav = 6.78*10**46
	print p
	## make sure R is consistent with Temp there!
	#qIR = (1./nu0)**(gam)
	#from scipy import special as spc
	#R = ma.sqrt(  Lav / (  4.* ma.pi * 8. * ma.pi  * qIR * h/c/c * (kb/h)**(4.+gam) * spc.gamma(4+gam) * spc.zetac(4.+gam) * Td**(4.+gam) ) )
	
	if (gam < 0):
		chi2 = np.inf
	else:
		#Rprint = sqtfR/pc2cm
		#print Rprint
		pref = np.ones(len(nu))
		for i in range(len(nu)):
			pref[i] = min(1., (nu[i]/nu0)**(gam))
		chi = (y - pref*Bv(nu, Td)* (sqtfR/Dst)**2 )/ dy
		chi2 = sum(chi*chi)

	print chi2
	return chi2








def Fsrc_Err2(p, t, y, dy, Args):
	print "EVAL", p
	#Lfac, bets, phs, incl = p
	Lfac, bets, phs, incl, Pday = p


	#FVbndRel, Lav, bets, Ombn, Dst = Args
	FVbndRel, Lav, Dst = Args

	if (bets<=0 or bets >=1.0):
		return inf
	else:
	#incl = np.arccos(0.067/bets)
		alphnu =1.1
		Ombn =	2.*ma.pi/(Pday*24.*3600.) * (1.+0.2784)
		chi = (y - -2.5*np.log10(Fsrc_Dop((t*3600.*24.-phs*2.*ma.pi/Ombn), Dst, ma.pi/2., 0.0, Lfac*Lav, bets, incl, Ombn, alphnu)/FVbndRel) )/ dy
		return sum(chi*chi)




def Fsrc_ISO_Err2(p, t, y, dy, Args):
	print "EVAL", p
	#Lfac, bets, phs, incl = p
	Lfac, Amp, phs, Pday = p

	#FVbndRel, Lav, bets, Ombn, Dst = Args
	FVbndRel, Lav, Dst = Args

	if (Amp<=0):
		return inf
	else:
	#incl = np.arccos(0.067/bets)
		Ombn =	2.*ma.pi/(Pday*24.*3600.) * (1.+0.2784)
		t0 = phs #* 2.*ma.pi/Ombn

		chi = (y - -2.5*np.log10(Fsrc_Iso((t*3600.*24.-t0), Dst, Lfac*Lav, Amp, Ombn, t0)/FVbndRel) )/ dy
		return sum(chi*chi)


##Optically Thick Torus Shell - Doppler Model
def magPoint_OpThick_TorShell(params, t, THEargs, RHStable, Ttable):

	sinJJ, cosTT, Rin = params

	n0 = 1000000.00000001
	Rin = Rin * 2.73213149e+18


	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
	
	FRel, numin, numax, Dist, Lav, Ombn, alph, pp, Rout,  aeff, nu0, nne, beta = THEargs
	IncFit = np.arccos(0.067/beta)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

	return -2.5*np.log10(F_ShTorOptThick_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)





def magPoint_OpThin_TorShell_dustP(params, t, THEargs, RHStable, Ttable):

	sinJJ, cosTT, Rin, nne, nu0 = params

	n0 = 1.0  #this shoudlnt matter opt thin is to IR, and is assumed in calcualtion method
	Rin = Rin * 2.73213149e+18

	aeff = (c/nu0)/(2.*ma.pi)

	if (nne <0.0):
		return np.inf

	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
	
	FRel, numin, numax, Dist, Lav, Ombn, alph, pp, Rout, beta = THEargs
	IncFit = np.arccos(0.067/beta)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

	return -2.5*np.log10(F_ShTorOptThin_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)



def magPoint_OpThin_TorShell(params, t, THEargs, RHStable, Ttable):

	#sinJJ, cosTT, Rin, alph = params
	sinJJ, cosTT, Rin = params

	alph = -2.0
	n0 = 1.0  #this shoudlnt matter opt thin is to IR, and is assumed in calcualtion method
	Rin = Rin * 2.73213149e+18




	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
	
	FRel, numin, numax, Dist, Lav, Ombn,     pp, Rout,  nu0, nne, beta = THEargs
	IncFit = np.arccos(0.067/beta)
	aeff = (c/nu0)/(2.*ma.pi)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

	return -2.5*np.log10(F_ShTorOptThin_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)



def magPoint_OpThin_TorShell_W1(params, t, THEargs, RHStable, Ttable):

	#sinJJ, cosTT, Rin, alph = params
	sinJJ, cosTT, RW1, RW2 = params

	alph = -2.0
	n0 = 1.0  #this shoudlnt matter opt thin is to IR, and is assumed in calcualtion method
	Rin = RW1 * pc2cm# 2.73213149e+18




	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
	
	FRel, numin, numax, Dist, Lav, Ombn,     pp, Rout,  nu0, nne, beta = THEargs
	IncFit = np.arccos(0.067/beta)
	aeff = (c/nu0)/(2.*ma.pi)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

	return -2.5*np.log10(F_ShTorOptThin_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)



def magPoint_OpThin_TorShell_W2(params, t, THEargs, RHStable, Ttable):

	#sinJJ, cosTT, Rin, alph = params
	sinJJ, cosTT, RW1, RW2 = params

	alph = -2.0
	n0 = 1.0  #this shoudlnt matter opt thin is to IR, and is assumed in calcualtion method
	Rin = RW2 * pc2cm# 2.73213149e+18




	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
	
	FRel, numin, numax, Dist, Lav, Ombn,     pp, Rout,  nu0, nne, beta = THEargs
	IncFit = np.arccos(0.067/beta)
	aeff = (c/nu0)/(2.*ma.pi)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

	return -2.5*np.log10(F_ShTorOptThin_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)





def OpThin_TorShell_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = (y1 - magPoint_OpThin_TorShell(p, t, THEargs1, RHStable, Ttable)) / dy1
	chi2 = (y2 - magPoint_OpThin_TorShell(p, t, THEargs2, RHStable, Ttable)) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2


def OpThin_TorShell_Err2_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = (y1 - magPoint_OpThin_TorShell_W1(p, t, THEargs1, RHStable, Ttable)) / dy1
	chi2 = (y2 - magPoint_OpThin_TorShell_W2(p, t, THEargs2, RHStable, Ttable)) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2




def magPoint_Sphere_dustP(params, t, THEargs, RHStable, Ttable):

	Rin, nne, nu0 = params

	n0 = 1.0  #this shouldn't matter opt thin is to IR, and is assumed in calcualtion method
	Rin = Rin * 2.73213149e+18

	aeff = (c/nu0)/(2.*ma.pi)

	if (nne <0.0):
		return np.inf

	t = t * 86400.
	JJ = 0.0#np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = 0.0
	
	FRel, numin, numax, Dist, Lav, Ombn, alph, pp, Rout, beta = THEargs
	IncFit = np.arccos(0.067/beta)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

	return -2.5*np.log10(F_ShTorOptThin_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)







def magPoint_OpThick_TorShell_dustP(params, t, THEargs, RHStable, Ttable):

	sinJJ, cosTT, Rin, nne, nu0 = params

	n0 = 1000000.00000001
	Rin = Rin * 2.73213149e+18

	aeff = (c/nu0)/(2.*ma.pi)

	if (nne <0.0):
		return np.inf

	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
	
	FRel, numin, numax, Dist, Lav, Ombn, alph, pp, Rout, beta = THEargs
	IncFit = np.arccos(0.067/beta)

	Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

	return -2.5*np.log10(F_ShTorOptThick_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)





	##Optically Thick Torus Shell - Doppler Model
def magPoint_OpThin_TorThick(params, t, THEargs, RHStable, Ttable):

	sinJJ, cosTT, Rin, n0, alph = params

	if (cosTT*cosTT > 1.0 or sinJJ*sinJJ > 1.0):
		return np.inf
	else:
		Rin = Rin * 2.73213149e+18
		n0 = n0*1.e-10
		t = t * 86400.
		JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
		thetT = np.arccos(cosTT)
		
		FRel, numin, numax, Dist, Lav, Ombn,     pp, Rout, aeff, nu0, nne, beta = THEargs
		IncFit = np.arccos(0.067/beta)

		Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

		return -2.5*np.log10(F_Thick_Dop_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)


          





def OpThick_TorShell_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = (y1 - magPoint_OpThick_TorShell(p, t, THEargs1, RHStable, Ttable)) / dy1
	chi2 = (y2 - magPoint_OpThick_TorShell(p, t, THEargs2, RHStable, Ttable)) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2


def OpThick_TorShell_dustP_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = (y1 - magPoint_OpThick_TorShell_dustP(p, t, THEargs1, RHStable, Ttable)) / dy1
	chi2 = (y2 - magPoint_OpThick_TorShell_dustP(p, t, THEargs2, RHStable, Ttable)) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2


def OpThin_TorShell_dustP_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = (y1 - magPoint_OpThin_TorShell_dustP(p, t, THEargs1, RHStable, Ttable)) / dy1
	chi2 = (y2 - magPoint_OpThin_TorShell_dustP(p, t, THEargs2, RHStable, Ttable)) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2



def Sphere_dustP_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = (y1 - magPoint_Sphere_dustP(p, t, THEargs1, RHStable, Ttable)) / dy1
	chi2 = (y2 - magPoint_Sphere_dustP(p, t, THEargs2, RHStable, Ttable)) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2



def DiffR_OpThick_TorShell_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	p1 = [p[0], p[1], p[2]]
	p2 = [p[0], p[1], p[3]]
	#p0 = [sinJ, cosT, Rin]
	chi1 = (y1 - magPoint_OpThick_TorShell(p1, t, THEargs1, RHStable, Ttable)) / dy1
	chi2 = (y2 - magPoint_OpThick_TorShell(p2, t, THEargs2, RHStable, Ttable)) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2



def OpThin_TorThick_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = (y1 - magPoint_OpThin_TorThick(p, t, THEargs1, RHStable, Ttable)) / dy1
	chi2 = (y2 - magPoint_OpThin_TorThick(p, t, THEargs2, RHStable, Ttable)) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2


	##Optically Thick Torus Shell - Doppler Model
def ISO_magPoint_OpThin_TorThick(params, t, THEargs, RHStable, Ttable):

	sinJJ, cosTT, Rin, n0, Amp = params


	if (cosTT*cosTT > 1.0 or sinJJ*sinJJ > 1.0):
		return np.inf
	else:
		Rin = Rin * 2.73213149e+18
		n0 = n0*1.e-10
		t = t * 86400.
		JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
		thetT = np.arccos(cosTT)
		
		FRel, numin, numax, Dist, Lav, Ombn, pp, aeff, nu0, nne, t0 = THEargs


		#Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
				#Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
		Aargs  = [Lav, Amp, Ombn, t0, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

		return -2.5*np.log10(F_Thick_Iso_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)

              

def ISO_OpThin_TorThick_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = ( y1 - ISO_magPoint_OpThin_TorThick(p, t, THEargs1, RHStable, Ttable) ) / dy1
	chi2 = ( y2 - ISO_magPoint_OpThin_TorThick(p, t, THEargs2, RHStable, Ttable) ) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2








def ISO_magPoint_OpThin_TorShell_W1(params, t, THEargs, RHStable, Ttable):

	#sinJJ, cosTT, Rin, Amp = params
	sinJJ, cosTT, RW1, RW2 = params

	Amp = 0.35

	n0 = 1.0 ## shouldnt matter
	if (cosTT*cosTT > 1.0 or sinJJ*sinJJ > 1.0):
		return np.inf
	else:
		Rin = RW1 * pc2cm#2.73213149e+18
		t = t * 86400.
		JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
		thetT = np.arccos(cosTT)
		
		FRel, numin, numax, Dist, Lav, Ombn, pp, aeff, nu0, nne, t0 = THEargs


		#Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
				#Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
		Aargs  = [Lav, Amp, Ombn, t0, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

		return -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)


def ISO_magPoint_OpThin_TorShell_W2(params, t, THEargs, RHStable, Ttable):

	#sinJJ, cosTT, Rin, Amp = params
	sinJJ, cosTT, RW1, RW2 = params

	Amp = 0.35

	n0 = 1.0 ## shouldnt matter
	if (cosTT*cosTT > 1.0 or sinJJ*sinJJ > 1.0):
		return np.inf
	else:
		Rin = RW2 * pc2cm#2.73213149e+18
		t = t * 86400.
		JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
		thetT = np.arccos(cosTT)
		
		FRel, numin, numax, Dist, Lav, Ombn, pp, aeff, nu0, nne, t0 = THEargs


		#Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
				#Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
		Aargs  = [Lav, Amp, Ombn, t0, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

		return -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)



def ISO_magPoint_OpThin_TorShell(params, t, THEargs, RHStable, Ttable):

	#sinJJ, cosTT, Rin, Amp = params
	sinJJ, cosTT, Rin = params

	Amp = 0.35

	n0 = 1.0 ## shouldnt matter
	if (cosTT*cosTT > 1.0 or sinJJ*sinJJ > 1.0):
		return np.inf
	else:
		Rin = Rin * 2.73213149e+18
		t = t * 86400.
		JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
		thetT = np.arccos(cosTT)
		
		FRel, numin, numax, Dist, Lav, Ombn, pp, aeff, nu0, nne, t0 = THEargs


		#Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
				#Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
		Aargs  = [Lav, Amp, Ombn, t0, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

		return -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)

              

def ISO_OpThin_TorShell_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = ( y1 - ISO_magPoint_OpThin_TorShell(p, t, THEargs1, RHStable, Ttable) ) / dy1
	chi2 = ( y2 - ISO_magPoint_OpThin_TorShell(p, t, THEargs2, RHStable, Ttable) ) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2



def ISO_OpThin_TorShell_Err2_TwoRs(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = ( y1 - ISO_magPoint_OpThin_TorShell_W1(p, t, THEargs1, RHStable, Ttable) ) / dy1
	chi2 = ( y2 - ISO_magPoint_OpThin_TorShell_W2(p, t, THEargs2, RHStable, Ttable) ) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2




	##Optically Thick Torus Shell - Doppler Model
def ISO_magPoint_OpThin_dustP_TorShell(params, t, THEargs, RHStable, Ttable):

	sinJJ, cosTT, Rin, Amp, nne, nu0 = params

	n0 = 1.0 ## shouldnt matter
	if (cosTT*cosTT > 1.0 or sinJJ*sinJJ > 1.0):
		return np.inf
	else:
		Rin = Rin * 2.73213149e+18
		t = t * 86400.
		JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
		thetT = np.arccos(cosTT)
		
		FRel, numin, numax, Dist, Lav, Ombn, pp, aeff, t0 = THEargs


		#Aargs  = [Lav, beta, IncFit, Ombn, alph, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]
				#Lavg, Amp, Ombin, t0, n0, Rd, p, thetT, JJ, aeff, nu0, nn = args
		Aargs  = [Lav, Amp, Ombn, t0, n0, Rin, pp, thetT, JJ, aeff, nu0, nne]

		return -2.5*np.log10(F_ShTorOptThin_Iso_QuadInt_PG(numin, numax, t, Dist, Aargs, RHStable, Ttable)/FRel)

              

def ISO_OpThin_TorShell_dustP_Err2(p, t, THEargs1, THEargs2, RHStable, Ttable, y1, dy1, y2, dy2):
	print "EVAL", p
	t1=time.clock()
	#p0 = [sinJ, cosT, Rin]
	chi1 = ( y1 - ISO_magPoint_OpThin_dustP_TorShell(p, t, THEargs1, RHStable, Ttable) ) / dy1
	chi2 = ( y2 - ISO_magPoint_OpThin_dustP_TorShell(p, t, THEargs2, RHStable, Ttable) ) / dy2
	sumChi2 = sum(chi1*chi1) + sum(chi2*chi2)
	print(sumChi2 )
	t2=time.clock()
	print(t2-t1)
	return sumChi2









