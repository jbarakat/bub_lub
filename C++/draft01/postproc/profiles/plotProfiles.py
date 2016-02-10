def load_src(name, fpath):
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def plotProfiles(STR, vL, epsL, CaL, xL, yL) :
	import matplotlib.pyplot as plt
	from matplotlib.pyplot import cm
	import numpy as np
	#from scipy import interpolate, interpolate
	from math import log10, floor, pi

	#epsL = round(epsL,-int(floor(log10(epsL)))+1)
	epsL = round(epsL,3)

	strvL  = str(vL)
	streL  = str(epsL)
	if (CaL >= 1) :
		CaL = int(CaL)
	strCL  = str(CaL)
	
	#yL = yL - np.mean(yL)
	
	plt.figure(figsize=(24,8))
	plt.plot(xL, yL, '-')
	
	# set plot preferences
	fs = 24
	plt.xlabel('$x/R$', fontsize=fs+4)

	if   (STR == 'SHAPE') :
		plt.ylabel('$y/R$',fontsize=fs+4)
		plt.ylim((-1,1))
	elif (STR == 'PRESS') :
		plt.ylabel('$pR / (\mu U)$',fontsize=fs+4)
	elif (STR == 'SHEAR') :
		plt.ylabel(r'$\tau R / (\mu U)$',fontsize=fs+4)
	elif (STR == 'TENS' ) :
		plt.ylabel(r'$\gamma / (\mu U)$',fontsize=fs+4)
	elif (STR == 'TRANS' ) :
		plt.ylabel(r'$q_s / (\mu U)$',fontsize=fs+4)
#	plt.xlim((-3.5,3.5))
	
	plt.xticks(fontsize=fs  )
	plt.yticks(fontsize=fs  )
	
	leg = plt.legend(('theory, ' + '$v$ = ' + strvL + ', $\epsilon$ = ' + streL + ', $\mathrm{Ca}$ = ' + strCL,),
	loc='best',fontsize=fs-4)
	leg.get_frame().set_linewidth(0.0)
	
	#plt.savefig(homedir + '/../plots/profiles/shape/v' + redvol + '/rmax' + maxrad + '/shape_v' + redvol + '_rmax' + maxrad + '_Ca' + capnum + '.pdf')
	#plt.close()
	plt.show()
