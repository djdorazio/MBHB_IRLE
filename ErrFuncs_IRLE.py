import numpy as np
import math as ma
#import numexpr as ne
import scipy as sc
from FluxFuncs_IRLE import *
import time




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





	##Optically Thick Torus Shell - Doppler Model
def magPoint_OpThin_TorThick(params, t, THEargs, RHStable, Ttable):

	sinJJ, cosTT, Rin, n0 = params

	Rin = Rin * 2.73213149e+18
	n0 = n0*1.e-10
	t = t * 86400.
	JJ = np.arcsin(sinJJ) ## CAREFUL WITH DOMAIN OF COS
	thetT = np.arccos(cosTT)
	
	FRel, numin, numax, Dist, Lav, Ombn, alph, Rout,  pp, aeff, nu0, nne, beta = THEargs
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





