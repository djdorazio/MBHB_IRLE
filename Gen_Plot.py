
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
import matplotlib.pyplot as plt

#import IR_LightEchoes_NewMeth as IRLE
from FluxFuncs_IRLE import *
from ErrFuncs_IRLE import *



def Plot_Shell_Thin_ISO(p, Nt,    tsrt, t_avg, t_MJD,    Lumsrt, W1_avsg, W2_avsg,   sigL, W1_sig, W2_sig):

	ttopt = np.linspace(tsrt[0]-100, t_MJD[len(t_MJD)-1]+100,       Nt)

	
	Fsrc_p_opt = [ 5.98144879e-02,   0.068,   6.55067929e-01,  -3.28334799e-04, 1.87091995e+03]


	
	OmFit = 2.*ma.pi/(Fsrc_p_opt[4]*24.*3600.)* (1.+0.2784)

	SET FROM DOP_PG
	opti = -2.5*np.log10(Fsrc_Dop_PG(ttopt*3600.*24, Dst, ma.pi/2., 0.0, Fsrc_p_opt[0]*Lav, Fsrc_p_opt[1], np.arccos(0.067/Fsrc_p_opt[1]), OmFit, 1.1)/FVbndRel)

	

	ttopt = (ttopt*(1.+zPG1302) - 50000)
	t_avg = (t_avg*(1.+zPG1302) - 50000)
	tsrt  = (tsrt*(1.+zPG1302) - 50000)
	t_MJD = (t_MJD*(1.+zPG1302) - 50000)






	plt.figure()
	plt.title("Doppler Source,  Torus Shell")
	plt.errorbar(tsrt, Lumsrt-3.0, yerr=sigL, linestyle="none", color = "blue", alpha=0.5) #alpha=0.1
	Fs = plt.plot(ttopt, opti-3.0, linestyle = '--', color='blue', linewidth=2)


	W1dat   = plt.errorbar(t_MJD, W1_mag, yerr=W1_sig, linestyle="none", color='orange', alpha=1., elinewidth=1.5)
	W2dat   = plt.errorbar(t_MJD, W2_mag+0.5, yerr=W2_sig, linestyle="none", color='red', alpha=1., elinewidth=1.5)

	W1av   = plt.errorbar(t_avg, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)
	W2av   = plt.errorbar(t_avg, W2_avg+0.5, yerr=W2_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)



	print "PLOT ISO Geo Thin OPT-THIN"
	W1shell = plt.plot(ttopt, magPoint_OpThin_TorShell(p, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
	W2shell = plt.plot(ttopt, magPoint_OpThin_TorShell(p , (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)


	plt.grid(b=True, which='both')
			

	plt.xlabel(r"$t$ [MJD]")
	plt.ylabel("mag")
	#plt.xlim(52000, 57500)
	plt.xlim(3000, max(ttopt))
	#plt.ylim(10.5, 11.5)
	plt.ylim(plt.ylim(10.5, 12.3)[::-1])

	#plt.show()
	plt.savefig("../emcee_data/"+Shell_File+"BestFit.png")
	plt.clf()





#(*SOME SYSTEM SPECIFIC CONSTANTS FOR TESTING*)
zPG1302 = 0.2784
Lav = 6.78*10**46 * 1.35
Dst = 1.4*10**9*pc2cm




## Wise band numbers
W1mx = numicron/2.8
W1mn = numicron/4.0
W2mx = numicron/3.9
W2mn = numicron/5.3

nuVbnd = c/(545*10**(-7))
FVbndRel = 3.636*10**(-20)*nuVbnd 
FW1Rel = 3.09540*10**(-20)*(W1mn + W1mx)/2
FW2Rel = 1.7187*10**(-20)*(W2mn + W2mx)/2





def Plot_Shell_Thin_ISO(p, Nt,    tsrt, t_avg, t_MJD,    Lumsrt, W1_avsg, W2_avsg,   sigL, W1_sig, W2_sig):


	ttopt = np.linspace(tsrt[0]-100, t_MJD[len(t_MJD)-1]+100,       Nt)

	
	Fsrc_ISO_p0 = [0.0597279747, 0.139181205, 0.688098413, 1871.99573]
	opti = -2.5*np.log10(Fsrc_Iso_PG(ttopt*3600.*24, Dst, Fsrc_p_opt[0]*Lav, Fsrc_ISO_p0[1], 1.0, 1.0)/FVbndRel)


	ttopt = (ttopt*(1.+zPG1302) - 50000)
	t_avg = (t_avg*(1.+zPG1302) - 50000)
	tsrt  = (tsrt*(1.+zPG1302) - 50000)
	t_MJD = (t_MJD*(1.+zPG1302) - 50000)






	plt.figure()
		plt.title("Isotropic Source, Torus Shell")
	plt.errorbar(tsrt, Lumsrt-3.0, yerr=sigL, linestyle="none", color = "blue", alpha=0.5) #alpha=0.1
	Fs = plt.plot(ttopt, opti-3.0, linestyle = '--', color='blue', linewidth=2)


	W1dat   = plt.errorbar(t_MJD, W1_mag, yerr=W1_sig, linestyle="none", color='orange', alpha=1., elinewidth=1.5)
	W2dat   = plt.errorbar(t_MJD, W2_mag+0.5, yerr=W2_sig, linestyle="none", color='red', alpha=1., elinewidth=1.5)

	W1av   = plt.errorbar(t_avg, W1_avg, yerr=W1_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)
	W2av   = plt.errorbar(t_avg, W2_avg+0.5, yerr=W2_avsg, linestyle="none", color='black', alpha=1., elinewidth=1.5)



	print "PLOT ISO Geo Thin OPT-THIN"
	W1shell = plt.plot(ttopt, ISO_magPoint_OpThin_TorShell(  p, (ttopt+50000)/(1.+zPG1302), W1args, RHS_table, T_table), linestyle = '--', color='orange', linewidth=2)
	W2shell = plt.plot(ttopt, ISO_magPoint_OpThin_TorShell(  p , (ttopt+50000)/(1.+zPG1302), W2args, RHS_table, T_table)+0.5, linestyle = '--', color='red', linewidth=2)


	plt.grid(b=True, which='both')
			

	plt.xlabel(r"$t$ [MJD]")
	plt.ylabel("mag")
	#plt.xlim(52000, 57500)
	plt.xlim(3000, max(ttopt))
	#plt.ylim(10.5, 11.5)
	plt.ylim(plt.ylim(10.5, 12.3)[::-1])

	#plt.show()
	plt.savefig("../emcee_data/"+Shell_File+"BestFit.png")
	plt.clf()