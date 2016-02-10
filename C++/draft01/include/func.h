/* FUNCTION EVALUATION
 *  Evaluate right-hand side vector f in the system of ODEs dy/dt = f(y,t).
 *
 * REFERENCES
 *  Secomb et al, J Fluid Mech (1986)
 *  
 * PARAMETERS
 *  ny					[input]		number of ODEs
 *  t						[input]		abscissas
 *  y						[input]		solution
 *  f						[input]		function
 *  r   = y[0]  [input]		radius
 *  psi = y[1]  [input]		tilt angle
 *  cs  = y[2]  [input]		meridional curvature
 *  qs  = y[3]  [input]		transverse shear tension
 *  p   = y[4]  [input]		pressure
 *  sig = y[5]  [input]		mean tension
 *  A   = y[6]  [input]		total surface area
 *  V   = y[7]  [input]		total volume
 *  Q   = y[8]  [input]		leakback flow rate
 *  S   = y[9]  [input]		total meridional arc length
 */

#ifndef FUNC_H
#define FUNC_H


/* HEADER FILES */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

using namespace std;

/* PROTOTYPES */

/* IMPLEMENTATIONS */
void func(int ny, double Ca,
          double t, double *y, double *f){
  int i;
  
  if (ny != 6)
    cout << "Error: ny should equal 6." << endl;
  
  // define variables
  double r   = y[0];
  double psi = y[1];
  double p   = y[2];
  double V   = y[3];
  double Q2  = y[4];
  double S   = y[5];
	
	// check bounds on radius
	if (r > 0.9999999999999)
		  r = 0.9999999999999;
	
	if (r <= 0)
		r = 1e-12;
	
	/*-------- RESULTS FROM LUBRICATION THEORY --------*/
	double r2  = r*r;
  double g;
  
	// calculate local axial pressure gradient, g = (dp/dx)/Ca
	if (r > 0.0001){
  	double log = gsl_sf_log(r  );

		g  = (8.0/(1.0 - r2));
		g *= (Q2 - (1.0 - r2));
		g /= (1.0 - 3.0*r2 - 4.0*r2*r2*log/(1.0 - r2));
	}
	else {
		g  = 8.0*(Q2 - 1.0); // limiting value as r --> 0
	}

	/*-------------------------------------------------*/
	
	if (r > 0.0001){ // check if far from end caps
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
  	double cosr = cos/r;
  	double sinr = sin/r;

  	// calculate function
  	f[0] = -sin;
  	f[1] = -p - cosr;
  	f[2] = Ca*g*cos;
  	f[3] = M_PI*r2*cos;
		f[4] = 0;
  	f[5] = 0;
	}
	else { // approximate cs, cphi as constant (spherical cap)
  	double cos = gsl_sf_cos(psi);
  	double sin = gsl_sf_sin(psi);
		
		f[0] = -sin;
		f[1] = -0.5*p;
		f[2] = Ca*g*cos; 
		f[3] = M_PI*r2*cos;
		f[4] = 0; 
		f[5] = 0;
	}

	for (i = 0; i < ny; i++){
		f[i] *= S;
	}
}



#endif
