import numpy as np
from numpy import *


### IR DATA
filename = "Short_Results.txt"

t_MJD = (np.genfromtxt(filename,usecols=23, comments="|"))/(1.+zPG1302) ## put in binary frame

W1_mag = np.genfromtxt(filename,usecols=7, comments="|")
W1_sig = np.genfromtxt(filename,usecols=8, comments="|")

W2_mag = np.genfromtxt(filename,usecols=10, comments="|")
W2_sig = np.genfromtxt(filename,usecols=11, comments="|")


W2_sig = np.genfromtxt(filename,usecols=11, comments="|")

##error flags
qual_frame  = np.genfromtxt(filename,usecols=19, comments="|")
saa_sep     = np.genfromtxt(filename,usecols=21, comments="|")
moon_masked = np.genfromtxt(filename,usecols=22, comments="|")

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