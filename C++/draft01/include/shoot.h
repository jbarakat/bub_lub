/* MULTIPLE SHOOTING METHOD
 *  Integrate ODEs using an IVP method over segments of the domain. Update
 *  solution using Newton's method by enforcing continuity and boundary
 *  conditions.
 *
 * REFERENCES
 *  Stoer & Bulirsch, pp. 557-561
 *  
 * PARAMETERS
 *  n			[input]		number of ODEs
 *  m			[input]		number of shooting points
 *	Ca		[input]		capillary number (parameter #1)
 *  area	[input]		total surface area (parameter #2)
 *  vlme	[input]		total volume (parameter #3)
 *  t			[input]		abscissa (ranges from 0 to 1)
 *  s			[input]		solution vector (size m*n)
 *  y			[input]		integral at interior points [size (m-1)*n]
 *  si		[input]		initial guess of solution
 *  sf		[output]	final solution after Newton iteration
 */

#ifndef SHOOT_H
#define SHOOT_H


/* HEADER FILES */
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <math.h>
#include "newton.h"

/* PROTOTYPES */
void mshoot(int, int, int, const int, const double, const double, double, double , double , double*, double*, double*, int&);
void shoot (int, int, int, double, double*, double*, double*);

/* IMPLEMENTATIONS */

/* Multiple shooting routine. Start with an initial guess si and
 * update guess using the Newton-Raphson method.
 */
void mshoot(int n, int m, int nrk, const int MAXITER, const double TOL, const double DTOL,
            double Ca, double area, double vlme,
						double *t, double *si, double *sf, int &flag){
	int    i, j, iter;
	int    k;
	int    mn  = m*n;
	int    m1n = (m-1)*n;
	int    info;
	double s[mn ], sp[mn ];
	double y[m1n], yp[m1n];
	double f[mn ], fp[mn ];
	double fnorm, f0norm, fpnorm, dfnorm;
	//const int MAXITER = 101;
	//const double TOL  = 1e-8;

	// parameters for Newton iteration
	double d[mn];
	double p = 1;
	
	/*------------- INITIALIZE -------------*/
	for (i = 0; i < mn; i++){
		s[i] = si[i];
	}
	iter   = 0;
	fnorm  = 1;
	f0norm = 0;
	dfnorm = 1;
	info   = 0;

	/*--------- TAKE MULTIPLE SHOTS --------*/
	while (iter < MAXITER && fnorm > TOL && dfnorm > DTOL && info == 0){
		shoot (n, m, nrk, Ca,             t, s, y   );
		newton(n, m, nrk, Ca, area, vlme, t, s, y, d, info);
		fzero (n, m, nrk, Ca, area, vlme, t, s, y, f);
		
		// calculate squared norm of f
		fnorm = 0;
		for (i = 0; i < mn; i++){
			fnorm += f[i]*f[i];
		}
		fnorm = sqrt(fnorm);

		/* use a finite search process to determine the smallest
		 * integer i that satisfies ||f(s+p*d)|| < ||f(s)||
		 * where p = 2^(-i). */
		for (i = 0; i < 5; i++){
			p = pow(2.0, -i);
			
			for (j = 0; j < mn; j++){
				sp[j] = s[j] - p*d[j];
			}

			shoot(n, m, nrk, Ca,             t, sp, yp    );
			fzero(n, m, nrk, Ca, area, vlme, t, sp, yp, fp);

			fpnorm = 0;
			for (j = 0; j < mn; j++){
				fpnorm += fp[j]*fp[j];
			}

			if (fpnorm < fnorm)
				break;
		}
		
		// update s
		for (i = 0; i < mn; i++){
			s[i] = sp[i];
		}
		
		cout << "iter = " << iter << ", fnorm = " << fnorm << ", p = " << p << endl;
		iter++;
	}

	// store solution
	for (i = 0; i < mn; i++){
		sf[i] = s[i];
	}

	if (fnorm < TOL*100)
		flag = 0;
	else
		flag = 1;

	dfnorm = fabs(fnorm - f0norm);
	f0norm = fnorm;
}

/* Shoot from both ends to the midpoint. Integrate the ODEs for y(t) using 
 * a guess s for the intermediate values of y at the grid points t.
 *
 * PARAMETERS:
 *  n = number of ODEs
 *  m = number of grid points
 *  t = an m-vector of grid points
 *  s = an mn-vector of the guess
 *  y = an (m-1)n-vector of the solution
 *    = [y(t1; t0, s0), y(t2; t1, s1), ..., y(tmh; tmh-1, smh-1), 
 *       y(tmh; tmh+1, smh+1), ..., y(tm-2; tm-1, sm-1)]
 */
void shoot(int n, int m, int nrk,
           double Ca, double *t, double *s, double *y){
	if (m % 2 == 0){
		cout << "Error: choose m to be odd." << endl;
		return;
	}
	
	int    i, j, i1;
	int    mh = (m+1)/2;
	double t0, t1;
	double s0[n], y1[n];

	// integrate forward in t
	for (i = 0; i < mh-1; i++){
		i1 = i + 1;

		// setup
		t0 = t[i ];
		t1 = t[i1];

		for (j = 0; j < n; j++){
			s0[j] = s[i *n + j];
		}

		// shoot from t0 to t1
		rk4(n, nrk, Ca, t0, t1, s0, y1);
		
		// store results
		for (j = 0; j < n; j++)
			y[i*n + j] = y1[j];
	}

	// integrate backward in t
	for (i = mh-1; i < m-1; i++){
		i1 = i + 1;

		// setup
		t0 = t[i1];
		t1 = t[i ];

		for (j = 0; j < n; j++){
			s0[j] = s[i1*n + j];
		}

		// shoot from t0 to t1
		rk4(n, nrk, Ca, t0, t1, s0, y1);
		
		// store results
		for (j = 0; j < n; j++)
			y[i*n + j] = y1[j];
	}
}

#endif
