def load_src(name, fpath) :
	import os, imp
	path = fpath if os.path.isabs(fpath) \
	else os.path.join(os.path.dirname(__file__), fpath)
	return imp.load_source(name, path)

def readProfiles() :
	import numpy as np
	#from scipy import interpolate, interpolate
	from math import log10, floor, pi, sqrt
	import os.path
	
	altdir = os.path.abspath("../")
	
	load_src("readLT", altdir + "/tools/readLT.py")
	from readLT import *
	
	###################################################
	#                     SETUP                       #
	###################################################
	
	# read parameters
	ifile = open('./params.in', 'r')
	print ifile
	for line in ifile :
		params = line.split()		# space delimited text
		if params[0] == 'PROFILE' :
			STR = params[2]
		if params[0] == 'REDVOL' :
			redvol = params[2]
		if params[0] == 'CONFIN' :
			confin = params[2]
		if params[0] == 'CAPNUM' :
			capnum = params[2]
	ifile.close()
		
	###################################################
	#               LUBRICATION THEORY                #
	###################################################
		
	xL = []
	yL = []

	vL   = 0
	epsL = 0
	CaL  = 0
	
	# load data file
	(Ca, eps, area, vlme, s, x, r, nx, nr, tx, tr, cs, cphi, qs, p, tau, gam) = readLT(altdir, redvol, confin, capnum)

#	lengthx = len(x)
#	mid = lengthx/2
#	print mid
#	print (1.0 - r[mid])/eps
#	print 1.0 - ((1.0 - r[mid])/eps - 1.0)*6.0/eps

	# reflect shape across symmetry axis
	xx  = np.concatenate((x, np.flipud( x)))
	rr  = np.concatenate((r, np.flipud(-r)))
	
	## interpolate to xi using splines
	#tck = interpolate.splrep(x, y, s=0)
	#yi = interpolate.splev(xi, tck, der=0)

	if   (STR == 'SHAPE') :
		xL = xx
		yL = rr
	elif (STR == 'PRESS') :
		xL = x
		yL = p
	elif (STR == 'SHEAR') :
		xL = x
		yL = tau
	#	for i in range(len(tau)) :
	#		yL[i] = tau[i]*Ca
	elif (STR == 'TENS' ) :
		xL = x
		yL = gam
	elif (STR == 'TRANS' ) :
		xL = x
		yL = qs
	
	# add a decimal point to the Ca string if need be
	if (capnum[0] == '0') :
		capnum = capnum[0] + '.' + capnum[1:]
	
	a  = sqrt(area/(4.0*pi))
	a3 = a*a*a
	v  = vlme/((4.0/3.0)*pi*a3)
	v  = round(v, 2)

	vL   = v
	epsL = eps
	CaL  = Ca

	return (STR, vL, epsL, CaL, xL, yL)
