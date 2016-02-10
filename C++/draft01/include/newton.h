/* NEWTON-RAPHSON METHOD
 *  Solve the linear system DF*d = F for the Newton step d.
 *
 * REFERENCES
 *  Stoer & Bulirsch, pp. 302-316, 559-561
 *  
 * PARAMETERS
 *  n		[input]			number of ODEs
 *  m		[input]			number of shooting points
 *  t		[input]			independent coordinate (ranges from 0 to 1)
 *  s		[input]			guess (an mn-vector)
 *  y		[input]			solution [an (m-1)n-vector]
 *  F		[output]		function to be zeroed
 *  DF	[output]		Jacobian
 */

#ifndef NEWTON_H
#define NEWTON_H


/* HEADER FILES */
#include <cstdlib>
#include <cstdio>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include "ode.h"
#include "bcond.h"

/* PROTOTYPES */
void newton(int, int, int, double, double, double,
						double *, double *, double *, double *, int&);
void fzero(int, int, int, double, double, double,
           double*, double*, double*, double*);
void jacob(int, int, int, double, double, double,
           double*, double*, double*, double*, double*);

/* IMPLEMENTATIONS */

/* Given last guess s and solution y, calculate the Newton step d. */
void newton(int n, int m, int nrk,
            double Ca, double  area, double vlme,
						double *t, double *s, double *y, double *d, int& info){
	int i, j, k, l;
	int i1, jj, kk;
	int ni;
	int mn   = m*n;
	int mnmn = mn*mn;
	int nn   = n*n;
	double F [mn], DF[mnmn];

	// size of Newton step (need to work on this for modified Newton's method)
	double p = 1;
	
	// generate F (function vector) and DF (Jacobian matrix)
	jacob(n, m, nrk, Ca, area, vlme, t, s, y, F, DF);
	
	/* solve the linear system using Gaussian elimination 
	 * (LAPACK's DGETRF, (DGEEQU,) and  DGETRS subroutines) */
	char   trans = 'N';				// form of the system of equations
	int    ldDF  =  mn;				// leading dimension of DF
	int    ldF   =  1;				// leading dimension of F
	int    nrhs  =  1;				// number of right-hand sides
	int    ipiv[mn];	 				// array of pivot indices

	// factorize
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, mn, mn, DF, ldDF, ipiv);
	if (info != 0){
		cout << "Error: dgetrf did not return successfully." << endl;
		return;
	}

	// solve
	info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, trans, mn, nrhs, DF, ldDF, ipiv, F, ldF);
	if (info != 0){
		cout << "Error: dgetrs did not return successfully." << endl;
		return;
	}
	
	// store solution
	for (i = 0; i < mn; i++){
		d[i] = F[i];
	}

//// IF CALCULATING THE JACOBIAN BY FINITE DIFFERENCES BECOMES TOO EXPENSIVE, 
//// IMPLEMENT BROYDEN'S UPDATE
//// -- note that Broyden shouldn't be used too many times in sequence
//// -- near the solution, Broyden's method becomes inaccurate
//// -- forward differences might be inaccurate near the solution as well..
////    (will need to use central difference)
}


/* Generate mn-vector F (function). */
void fzero(int n, int m, int nrk,
           double Ca, double area, double vlme,
           double *t, double *s, double *y, // takes t, s, and y as input
           double *F){
  if (m % 2 == 0){
    cout << "Error: choose m to be odd." << endl;
    return;
  }

  int i, j, k, l;
  int jj, kk, ll;
  int mh = (m+1)/2;
  int mn = m*n;
  int nn = n*n;
  int ni, ni1, ni2;
  double A[nn], B[nn], c[n], r[n];
  double s0,  s1;
  double t0,  t1;

  // calculate linear BCs
  bcond(n, area, vlme, A, B, c);
  for (i = 0; i < n; i++){
    r[i] = -c[i];
    for (j = 0; j < n; j++){
      k  = i*n + j;
      kk = (m-1)*n + j;
      s0 = s[j];
      s1 = s[kk];

      r[i] += A[k]*s0 + B[k]*s1;
    }
  }

  /*--------- CALCULATE FUNCTION VECTOR ---------*/
  // continuity conditions
  for (i = 0; i < mh-1; i++){
    for (j = 0; j < n; j++){
      k  = i*n + j;
      kk = (i+1)*n + j;

      F[k] = y[k] - s[kk];
    }
  }

  for (i = mh-1; i < m-1; i++){
    for (j = 0; j < n; j++){
      k  = i*n + j;

      F[k] = y[k] - s[k];
    }
  }

  // boundary conditions
  i = m-1;
  for (j = 0; j < n; j++){
    k = i*n + j;
    F[k] = r[j];
  }
}

void jacob(int n, int m, int nrk,
           double Ca, double area, double vlme,
           double *t, double *s, double *y,
           double *F, double *DF){
  if (m % 2 == 0){
    cout << "Error: choose m to be odd." << endl;
    return;
  }

  int i, j, k, l;
  int jj, kk, ll;
  int mh = (m+1)/2;
  int mn = m*n;
  int nn = n*n;
  int ni, ni1, ni2;
  double A[nn], B[nn], c[n], r[n], G[nn];
  double sp[n], yp[n];
  double s0, s1;
  double t0, t1;
  double Ds = 0.0001;
  double Dy;

  // calculate linear BCs
  bcond(n, area, vlme, A, B, c);
  for (i = 0; i < n; i++){
    r[i] = -c[i];
    for (j = 0; j < n; j++){
      k  = i*n + j;
      kk = (m-1)*n + j;
      s0 = s[j];
      s1 = s[kk];

      r[i] += A[k]*s0 + B[k]*s1;
    }
  }
 
  /*--------- CALCULATE FUNCTION VECTOR ---------*/
  // continuity conditions
  for (i = 0; i < mh-1; i++){
    for (j = 0; j < n; j++){
      k  = i*n + j;
      kk = (i+1)*n + j;

      F[k] = y[k] - s[kk];
    }
  }
 
  for (i = mh-1; i < m-1; i++){
    for (j = 0; j < n; j++){
      k  = i*n + j;

      F[k] = y[k] - s[k];
    }
  }
	
  // boundary conditions
  i = m-1;
  for (j = 0; j < n; j++){
    k = i*n + j;
    F[k] = r[j];
  }


  /*--------- CALCULATE JACOBIAN MATRIX ---------*/
  /*--------- USING FORWARD DIFFERENCES ---------*/
  // initialize
  for (i = 0; i < mn; i++){
    for (j = 0; j < mn; j++){
      DF[i*mn + j] = 0;
    }
  }

  // assign G and -I (first mh-1 blocks)
  for (i = 0; i < mh-1; i++){   // block of n rows
    ni  = i*n;
    ni1 = (i+1)*n;
    ni2 = (i+2)*n;
    t0  = t[i  ];
    t1  = t[i+1];

    // assign G to main diagonal
    for (k = 0; k < n; k++){    // columns: ds0 to ds(n-1)
      kk = k + ni;

      // perturb sp = s + Ds
      for (j = 0; j < n; j++){
        jj = j + ni;
        sp[j] = s[jj];
      }
      sp[k] += Ds;

//      // check radius
//      if (k == 0 && sp[k] > 0.99999999999)
//        sp[k] -= 2.0*Ds;

      // check tilt angle
      if (k == 1 && sp[k] > M_PI/2)
        sp[k] -= 2.0*Ds;

      // integrate to get yp
      rk4(n, nrk, Ca, t0, t1, sp, yp);

      /* calculate forward difference part
       * and update G */
      for (j = 0; j < n; j++){  // rows: dy0 to dy(n-1)
        jj = j + ni;
        Dy = yp[j] - y[jj];
      //  G[j*n  + k] = 0.5*Dy/Ds;
        G [j *n   + k ] = Dy/Ds;
        DF[jj*mn  + kk] = G[j*n + k];
      }

    //  // perturb sp = s - Ds
    //  sp[k] -= 2*Ds;

    //  // check radius
    //  if (k == 0 && sp[k] < 0)
    //    sp[k] += 2.0*Ds;

    //  // check tilt angle
    //  if (k == 1 && sp[k] < -M_PI/2)
    //    sp[k] += 2.0*Ds;

    //  // integrate to get yp
    //  rk4(n, nrk, Ca, t0, t1, sp, yp);

    //  /* calculate backward difference part
    //   * and update G and DF */
    //  for (j = 0; j < n; j++){  // rows: dy0 to dy(n-1)
    //    jj = j + ni;
    //    Dy = y[jj] - yp[j];
    //  
    //    G [j *n  + k ] += 0.5*Dy/Ds;
    //    DF[jj*mn + kk]  = G[j*n + k];
    //  }
    }

    // assign -I to upper diagonal
    for (j = ni; j < ni1; j++){
      for (k = ni1; k < ni2; k++){
        if (j == k-n){
          DF[j*mn + k] = -1;
        }
      }
    }
  }

  // assign G and -I (second mh-1 blocks)
  for (i = mh-1; i < m-1; i++){ // block of n rows
    ni  = i*n;
    ni1 = (i+1)*n;
    ni2 = (i+2)*n;
    t0  = t[i+1];
    t1  = t[i  ];

    // assign -I to main diagonal
    for (j = ni; j < ni1; j++){
      for (k = ni; k < ni1; k++){
        if (j == k){
          DF[j*mn + k] = -1;
        }
      }
    }

    // assign G to upper diagonal
    for (k = 0; k < n; k++){    // columns: ds0 to dsn
      kk = k + ni1;

      // perturb sp = s + Ds
      for (j = 0; j < n; j++){
        jj = j + ni1;
        sp[j] = s[jj];
      }
      sp[k] += Ds;

      // check radius
//      if (k == 0 && sp[k] > 0.999999999)
//        sp[k] -= 2.0*Ds;

      // check tilt angle
      if (k == 1 && sp[k] > M_PI/2)
        sp[k] -= 2.0*Ds;

      // integrate to get yp
      rk4(n, nrk, Ca, t0, t1, sp, yp);

      /* calculate forward difference part
       * and update G */
      for (j = 0; j < n; j++){  // rows: dy0 to dyn
        jj = j + ni;
        Dy = yp[j] - y[jj];

    //    G[j*n + k] = 0.5*Dy/Ds;
        G [j *n  + k ] = Dy/Ds;
        DF[jj*mn + kk] = G[j*n + k];
      }

    //  // perturb sp = s - Ds
    //  sp[k] -= 2*Ds;

    //  // check radius
    //  if (k == 0 && sp[k] < 0)
    //    sp[k] += 2.0*Ds;

    //  // check tilt angle
    //  if (k == 1 && sp[k] < -M_PI/2)
    //    sp[k] += 2.0*Ds;

    //  // integrate to get yp
    //  rk4(n, nrk, Ca, t0, t1, sp, yp);

    //  /* calculate backward difference part
    //   * and update G and DF*/
    //  for (j = 0; j < n; j++){  // rows: dy0 to dyn
    //    jj = j + ni;
    //    Dy = y[jj] - yp[j];
    //  
    //    G [j *n  + k ] += 0.5*Dy/Ds;
    //    DF[jj*mn + kk]  = G[j*n + k];
    //  }
    }
  }

  // assign A and B to bottom left and right corners
  i   = m-1;
  ni  = i*n;
  for (j = 0; j < n; j++){
    jj = j + ni;
    for (k = 0; k < n; k++){
      kk = k + ni;
      DF[jj*mn + k ] = A[j*n + k];
      DF[jj*mn + kk] = B[j*n + k];
    }
  }
}

#endif
