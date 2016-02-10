/* WRITE FILE
 *  Write .dat file containing the abscissa and solution vectors.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  n			[input]		number of ODEs
 *  m			[input]		number of shooting points
 *  v			[input]		reduced volume (x 100)
 *  conf	[input]		confinement parameter scaled by critical value (x 100)
 *  t			[output]	abscissas (meridional arc length scaled to the [0,1] axis)
 *  r			[output]	cylindrical radius
 *  psi		[output]	tilt angle
 *  cs		[output]	meridional curvature (points outward from closed contour)
 *  qs		[output]	meridional component of transverse shear tension
 *  p			[output]	pressure difference between exterior and interior
 *  sig		[output]	mean tension
 *  A			[output]	surface area
 *  V			[output]	volume
 *	Q2		[output]	leakback flow rate
 *  S			[output]	total meridional arc length (half-space)
 */

#ifndef WRITE_H
#define WRITE_H

/* HEADER FILES */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */

/* Write .dat file. */
void writeSoln(int n, int m, int v, double Ca, 
               double *t, double *s, string path){
	// error flag
	if (n != 6){
		cout << "Error: support only for n = 6." << endl;
		return;
	}

	// declare variables
	int i;
	int width = 12;
	ostringstream ssRedVol, ssConfin, ssCapNum;
	ostringstream header;
	ofstream file;
	string dir, fn, line;

	// set transverse shear tension to zero everywhere
	vector<double> QS(m, 0.0);

	// calculate meridional curvature using Laplace relation
	vector<double> CS(m);
	double rr, ppsi, pp, ssig;
	for (i = 0; i < m; i++){
		rr   = s[i*n + 0];
		ppsi = s[i*n + 1];
		pp   = s[i*n + 2];
		if (rr > 0.002)
			CS[i] = pp + gsl_sf_cos(ppsi)/rr;
		else
			CS[i] = 0.5*pp;
	}
	
	// calculate area using finite differences
	vector<double> AREA(m, 0.0);
	int i1;
	double r1, S1, da;
	for (i = 0; i < m-1; i++){
		if (i == 0)
			continue;
		else
			i1 = i + 1;
			r1 = s[i1*n + 0];
			S1 = s[i1*n + 5];
			da = 2.0*M_PI*r1*S1*(t[i1] - t[i]);
			AREA[i1] = AREA[i] + da;
	}

	// convert integers to strings
	ssRedVol << v   ;

	// convert double to string
	if (Ca >= 0.000000001 && Ca < 0.00000001)
		ssCapNum << fixed << setprecision(10) << Ca;
	if (Ca >= 0.00000001 && Ca < 0.0000001)
		ssCapNum << fixed << setprecision(9) << Ca;
	if (Ca >= 0.0000001 && Ca < 0.000001)
		ssCapNum << fixed << setprecision(8) << Ca;
	if (Ca >= 0.000001 && Ca < 0.00001)
		ssCapNum << fixed << setprecision(7) << Ca;
	if (Ca >= 0.00001 && Ca < 0.0001)
		ssCapNum << fixed << setprecision(6) << Ca;
	if (Ca >= 0.0001 && Ca < 0.001)
		ssCapNum << fixed << setprecision(5) << Ca;
	if (Ca >= 0.001 && Ca < 0.01)
		ssCapNum << fixed << setprecision(4) << Ca;
	if (Ca >= 0.01 && Ca < 0.1)
		ssCapNum << fixed << setprecision(3) << Ca;
	if (Ca >= 0.1 && Ca < 1.0)
		ssCapNum << fixed << setprecision(2) << Ca;
	if (Ca >= 1.0)
		ssCapNum << fixed << setprecision(0) << Ca;
	
	// set directory
	dir = path + "/v"    + ssRedVol.str(); 
	
	// set filename
	fn = "/sln_v" + ssRedVol.str()
		 + "_Ca"    + ssCapNum.str();
	fn.erase(remove(fn.begin(), fn.end(), '.'), fn.end());
	fn = dir + fn + ".dat";

	// write to file
	file.open(fn.c_str());
	header << setw(width)   << "t" 
	       << setw(width)   << "r" 
				 << setw(width)   << "psi"
				 << setw(width)   << "cs"
				 << setw(width)   << "qs"
	       << setw(width+4) << "p" 
				 << setw(width+4) << "sig"
				 << setw(width)   << "A"
				 << setw(width)   << "V"
				 << setw(width)   << "Q2"
				 << setw(width)   << "S";
	file << header.str() << endl;
	for (i = 0; i < m; i++){
		ostringstream tt, r, psi, cs, qs, p, 
		                  sig, A, V , Q2, S;
		tt  << fixed << setw(width)   << t[i];
		r   << fixed << setw(width)   << s[i*n + 0];
		psi << fixed << setw(width)   << s[i*n + 1];
		cs  << fixed << setw(width)   << CS[i];
		qs  << fixed << setw(width)   << QS[i];
		p   << fixed << setw(width+4) << s[i*n + 2];
		sig << fixed << setw(width+4) << 1.0;
		A   << fixed << setw(width)   << AREA[i];
		V   << fixed << setw(width)   << s[i*n + 3];
		Q2  << fixed << setw(width)   << s[i*n + 4];
		S   << fixed << setw(width)   << s[i*n + 5];
		
		line = tt.str() + r.str() + psi.str() + cs.str() + qs.str()
		     + p.str() + sig.str() + A.str() + V .str() + Q2.str() + S.str();

		file << line << endl;
	}
	file.close();

	cout << "Solution written to " << fn << "." << endl;
}


#endif
