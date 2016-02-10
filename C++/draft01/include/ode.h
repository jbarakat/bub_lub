/* ORDINARY DIFFERENTIAL EQUATIONS
 *  Integration schemes for a system of ODEs.
 *
 * REFERENCES
 *  Moin, Cambridge University Press (2010), pp. 64-70
 *  
 * PARAMETERS
 *  ny		[input]		number of elements in y
 *  nstep	[input]		number of integration steps
 *  p			[input]		parameter
 *  t0		[input]		lower boundary of domain
 *  t1		[input]		upper boundary of domain
 *  y0		[input]		solution vector at t = t0 (initial condition)
 *  y1		[output]	solution vector at t = t1
 */

#ifndef ODE_H
#define ODE_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include "func.h"

/* PROTOTYPES */
void rk4(int, int, double, double, double, double *, double *);
void eulF(int, int, double, double, double, double *, double *);

/* IMPLEMENTATIONS */

/* Given t0 and y0 (with parameter p), integrate to y1 at t1
 * using a fourth-order Runge-Kutta scheme. */
void rk4(int ny, int nstep, double p,
         double t0, double t1, double *y0, double *y1){
  int istep, iy, j;
  double t;
  double dt = (t1 - t0)/nstep;
  double *y00, *y01, *y02;
  double *f;

  y00 = (double*) calloc(ny , sizeof(double));
  y01 = (double*) calloc(ny , sizeof(double));
  y02 = (double*) calloc(ny , sizeof(double));
  f   = (double*) calloc(ny , sizeof(double));
  
  // initialize
  t = t0;
  for (iy = 0; iy < ny; iy++)
    y02[iy] = y0[iy];

  for (istep = 0; istep < nstep; istep++){
//		cout << istep << endl;

    // setup
    t += istep*dt;
    for (iy = 0; iy < ny; iy++)
      y00[iy] = y02[iy];

    // first step
    for (iy = 0; iy < ny; iy++)
      y01[iy] = y00[iy];
    func(ny, p, t, y01, f);
    for (iy = 0; iy < ny; iy++)
      y02[iy] += dt*(1.0/6.0)*f[iy];

    // second and third steps
    for (j = 0; j < 2; j++){
      for (iy = 0; iy < ny; iy++)
        y01[iy] = y00[iy] + dt*0.5*f[iy];
      func(ny, p, t, y01, f);
      for (iy = 0; iy < ny; iy++)
        y02[iy] += dt*(1.0/3.0)*f[iy];
    }

    // fourth step
    for (iy = 0; iy < ny; iy++)
      y01[iy] = y00[iy] + dt*f[iy];
    func(ny, p, t, y01, f);
    for (iy = 0; iy < ny; iy++)
      y02[iy] += dt*(1.0/6.0)*f[iy];
  }

  // store results
  for (iy = 0; iy < ny; iy++)
    y1[iy] = y02[iy];
  
  free(y00);
  free(y01);
  free(y02);
  free(f );
}

/* Given t0 and y0 (with parameter p), integrate to y1 at t1
 * using a forward Euler scheme. */
void eulF(int ny, int nstep, double p,
         double t0, double t1, double *y0, double *y1){
  int istep, iy, j;
  double t;
  double dt = (t1 - t0)/nstep;
	double *y00;
  double *f;
  
  y00 = (double*) malloc(ny * sizeof(double));
  f   = (double*) malloc(ny * sizeof(double));

	// initialize
	for (iy = 0; iy < ny; iy++)
		y00[iy] = y0[iy];
	
	for (istep = 0; istep < nstep; istep++){
//		cout << istep << endl;

		// setup
		t += istep*dt;

		// take step
		func(ny, p, t, y00, f);
		for (iy = 0; iy < ny; iy++)
			y00[iy] += dt*f[iy];
	}

	// store results
	for (iy = 0; iy < ny; iy++)
		y1[iy] = y00[iy];

	free(y00);
	free(f) ;
}

#endif
