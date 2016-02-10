def convex_hull(xy) :
	import matplotlib.pyplot as plt
	import numpy as np
	
	# sort data by x-coordinate
	xy = sorted(xy, key=lambda x:x[0])
	xy = np.column_stack(xy)
	x = xy[0]
	y = xy[1]
	n = len(x)
	
	# initialize arrays
	Ux = []
	Uy = []
	Lx = []
	Ly = []
	
	# append to Ux, Uy (upper hull)
	Ux.append(x[0	])	
	Ux.append(x[1	])	
	Uy.append(y[0	])	
	Uy.append(y[1	])	
	nU = len(Ux)
	for i in range(n-2) :
		j = i + 2 
		Ux.append(x[j])
		Uy.append(y[j])
		nU = len(Ux)
		
		# check if the last three points make a right turn
		# (if so, delete the middle point)
		while (nU > 2) :
			px = Ux[nU-3] - Ux[nU-1]
			py = Uy[nU-3] - Uy[nU-1]
			qx = Ux[nU-2] - Ux[nU-1]
			qy = Uy[nU-2] - Uy[nU-1]
			c	= px*qy - py*qx
			if (c >= 0) :
				Ux.pop(nU-2)
				Uy.pop(nU-2)
				nU = len(Ux)
			else :
		 		break
	
	# append to Lx, Ly (lower hull)
	Lx.append(x[n-1])
	Lx.append(x[n-2])
	Ly.append(y[n-1])
	Ly.append(y[n-2])
	nL = len(Lx)
	for i in range(n-2) :
		j = i + 3
		k = n - j
		Lx.append(x[k])
		Ly.append(y[k])
		nL = len(Lx)
	
		# check if the last three points make a right turn
		# (if so, delete the middle point)
		while (nL > 2) :
			px = Lx[nL-1] - Lx[nL-3]
			py = Ly[nL-1] - Ly[nL-3]
			qx = Lx[nL-2] - Lx[nL-3]
			qy = Ly[nL-2] - Ly[nL-3]
			c	= px*qy - py*qx
			if (c <= 0) :
				Lx.pop(nL-2)
				Ly.pop(nL-2)
				nL = len(Lx)
			else :
				break
	
	# remove the first and last elements from Lx, Ly
	# (to avoid duplication of the points from Ux, Uy)
#	Lx.pop(nL-1)
#	Ly.pop(nL-1)
#	Lx.pop(0)
#	Ly.pop(0)
	nL = len(Lx)
	
	# concatenate arrays and return result
	Hx = np.concatenate((Ux,Lx))
	Hy = np.concatenate((Uy,Ly))	
	H	= np.column_stack([Hx,Hy])
	return H
