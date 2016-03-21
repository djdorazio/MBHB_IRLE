import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt

import IR_LightEchoes_fast as IRLE
from IR_LightEchoes_fast import *


import numpy as np
import math as ma
import numexpr as ne
import scipy as sc
from scipy.optimize import brentq #fmin

import scipy.integrate as intg



def QvBv(nu, T, nu0, nn):
	qv = (nu/nu0)**(nn)
	if (type(qv) is float):
		if (qv>1.0):
			qv=1.0
	else:
		ii = np.where(qv > 1.0)[0]
		qv[ii] = 1.0
	return 2.*h*nu*nu*nu/(c*c)*1./(ma.exp(h*nu/(kb*T)) - 1.) * qv

## equation to minimize to solve for Temp
def T_RHS(Td, nu0, nn):
	RHS = 4.*ma.pi*(intg.quad(QvBv  ,0., nu0 , args=(Td, nu0, nn) )[0] + intg.quad(Bv  ,nu0 ,numicron*1000., args=(Td) )[0])#, epsabs=myabs, epsrel=myrel, limit=reclim, limlst = limlst, maxp1=maxp1, full_output=fo  )[0]
	return np.log10(RHS)



#must be?
nne = 1.
nu0 = numicron*0.2

NT = 10000
RHS_table = np.zeros(NT)
T_table = np.linspace(200., 1800., NT)
for i in range(NT):
	RHS_table[i] = T_RHS(T_table[i], nu0, nne)



