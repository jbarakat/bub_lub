/* SPHEROID
 *  Auxiliary functions to generate prolate spheroids of given area and volume
 *  and, conversely, calculate the area and volume of a spheroid of given
 *  major and minor radii.
 *
 *  NOTE (8/31/2015): This class is not complete. For this class to be truly 
 *  useful, should also write a function to calculate the curvature of the
 *  spheroid along the meridional contour.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  a			[input/output]		minor radius
 *	b			[input/output]		major radius
 *  area	[input/output]		total surface area
 *  vlme	[input/output]		total volume
 */

#ifndef SPHEROID_H
#define SPHEROID_H

/* HEADER FILES */
#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

/* PROTOTYPES */
int    proAxes(double, double, double, double, double &, double &);
double proArea(double, double);
double proVlme(double, double);

/* IMPLEMENTATIONS */

int proAxes(double area, double vlme, double a0, double b0, double &a, double &b){
	int    i, j, iter;

	double s[2], f[2], Df[4], Ds[2];
	double detDf, invDf[4];
	double fnorm, fpnorm;
	double A, V;
	double dAda, dAdb;
	double dVda, dVdb;
	double a2, b2;
	double e2, d2;
	double e , d ;
	double asine ;
	double p, sp[2], fp[2];

	const double PI2 = 2.0*M_PI;
	const int MAXITER = 1000;
	const double TOL  = 1e-16;
	
	// initial guess
	a    = a0;
	b    = b0;
	s[0] = a;
	s[1] = b;
	
	iter  = 0;
	fnorm = 1;
	while (iter < MAXITER && fnorm > TOL){
		if (a > b){
			cout << "Error: a > b. Pick another initial guess." << endl;
			return (1);
		}
		if (a < 0){
			cout << "Error: a < 0. Pick another initial guess." << endl;
			return (1);
		}
		if (b < 0){
			cout << "Error: b < 0. Pick another initial guess." << endl;
			return (1);
		}

		a2    = a*a;
		b2    = b*b;
		e2    = 1.0 - a2/b2;
		e     = sqrt(e2);
		d2    = 1.0 - e2;
		d     = sqrt(d2);
		asine = gsl_complex_abs(gsl_complex_arcsin_real(e));

		A  = proArea(a, b);
		V  = proVlme(a, b);
		
		f[0] = A - area; f[1] = V - vlme;

		fnorm = f[0]*f[0] + f[1]*f[1];
		
		/* calculate Jacobian
		 *  Df = / dA/da  dA/db \
		 *       \ dV/da  dV/db / */
		dAda  = -a/e2;
		dAda +=  a2*asine/(b*e2*e);
		dAda +=  b*asine/e;
		dAda +=  2.0*a;
		dAda *=  PI2;
	
		dAdb  =  a2/(b*e2);
		dAdb +=  a*asine/e;
		dAdb -=  a*a2*asine/(b2*e2*e);
		dAdb *=  PI2;
	
		dVda  =  4.0*PI2*a*b/3.0;
		dVdb  =  2.0*PI2*a*a/3.0;
	
		Df[0] = dAda; Df[1] = dAdb;
		Df[2] = dVda; Df[3] = dVdb;
	
		// invert Df and calculate Newton step Ds
		invDf[0] =  Df[3]; invDf[1] = -Df[1];
		invDf[2] = -Df[2]; invDf[3] =  Df[0];
		detDf = Df[0]*Df[3] - Df[1]*Df[2];
	
		for (i = 0; i < 4; i++){
			invDf[i] /= detDf;
		}
	
		for (i = 0; i < 2; i++){
			Ds[i] = 0;
			for (j = 0; j < 2; j++){
				Ds[i] += -invDf[i*2 + j]*f[j];
			}
		}

		// line search and backtracking
		p = 1;
		for (i = 0; i < 4; i++){
			p = pow(2.0, -i);

			sp[0] = s[0] + p*Ds[0];
			sp[1] = s[1] + p*Ds[1];

			A = proArea(sp[0], sp[1]);
			V = proVlme(sp[0], sp[1]);

			fp[0] = A - area;
			fp[1] = V - vlme;

			fpnorm = f[0]*f[0] + f[1]*f[1];

			if (fpnorm < fnorm)
				break;
		}

	//	cout << Ds[0] << " " << Ds[1] << endl;
		
		// update next guess
		for (i = 0; i < 2; i++){
			s[i] += p*Ds[i];
		}
	
		a = s[0]; b = s[1];

		cout << "iter = " << iter << ", p = " << p << ", fnorm = " << fnorm << endl;
		iter++;
	}

	return(0);
}

double proArea(double a, double b){ // b > a
	double area;
	double a2    = a*a;
	double b2    = b*b;
	double e2    = 1.0 - a2/b2;
	double e     = sqrt(e2);
	double asine = gsl_complex_abs(gsl_complex_arcsin_real(e));

	area  = 1 + b/(a*e)*asine;
	area *= 2*M_PI*a2;

	return(area);
}

double proVlme(double a, double b){ // b > a
	double vlme;
	double a2 = a*a;

	vlme  = a2*b;
	vlme *= 4.0*M_PI/3.0;

	return(vlme);
}


// UNFINISHED FUNCTION
void proShape(int m, double a, double b, 
						  double &S, double *t, double *r, double *cs, double *psi){ // b > a
	int    i, i1;
	double s[m], x[m], thet[m];
	double a2 = a*a;
	double b2 = b*b;
	double x2, r2;
	double dx, dx2;
	double dr, dr2;
	double ds, ds2;
	double s0 = 0;
	
	// get axial length
	x[0] = -b;
	dx   = 2*b/(m-1);
	for (i = 0; i < m; i++){
		i1 = i + 1;
		x[i1] = x[i] + dx;
	}

	// get radius
	for (i = 0; i < m; i++){
		x2   = x[i]*x[i];
		r2   = a2*(1.0 - x2/b2);
		r[i] = sqrt(r2);
	}

	// get arclength using forward difference method
	s[0] = s0;
	for (i = 0; i < m-1; i++){
		i1    = i + 1;
		dx    = x[i1] - x[i]; dx2 = dx*dx;
		dr    = r[i1] - r[i]; dr2 = dr*dr;
		ds2   = dx2 + dr2;
		ds    = sqrt(ds2);
		s[i1] = s[i] + ds;
	}

	S = s[m-1];

	// get t parameter
	for (i = 0; i < m; i++){
		t[i] = s[i]/S;
	}
	
	// get polar angle
	for (i = 0; i < m; i++){
		thet[m] = gsl_complex_abs(gsl_complex_arccos_real(x[i]/b));
	}

	// UNFINISHED (8/31/2015)




}

#endif
