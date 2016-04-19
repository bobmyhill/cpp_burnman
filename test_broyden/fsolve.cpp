/*  broyden.cpp -- Broyden's method for solving a nonlinear system of
 *              equations.
 * Modified from the code downloaded from http://www.crbond.com
 *  (C) 2001, C. Bond. All rights reserved.
 */
#include <cmath>
#include "fsolve.hpp"

  /* This version preserves the matrix 'a' and the vector 'b'. */
int fsolve::gelimd2(double **a,double *b,double *x, int n)
{
  double tmp,pvt,*t,**aa,*bb;
  int i,j,k,retval;
  
  retval = 0;
  aa = new double *[n];
  bb = new double [n];
  for (i=0;i<n;i++) {
    aa[i] = new double [n];
    bb[i] = b[i];
    for (j=0;j<n;j++) {
      aa[i][j] = a[i][j];
    }
  }
  for (i=0;i<n;i++) {             // outer loop on rows
    pvt = aa[i][i];              // get pivot value
    if (!pvt) {
      for (j=i+1;j<n;j++) {
	if((pvt = a[j][i]) != 0.0) break;
      }
      if (!pvt) {
	retval = 1;
	goto _100;     // pull the plug!
      }
      t=aa[j];                 // swap matrix rows...
      aa[j]=aa[i];
      aa[i]=t;
      tmp=bb[j];               // ...and result vector
      bb[j]=bb[i];
      bb[i]=tmp;        
    }
    // (virtual) Gaussian elimination of column
    for (k=i+1;k<n;k++) {       // alt: for (k=n-1;k>i;k--)
      tmp = aa[k][i]/pvt;
      for (j=i+1;j<n;j++) {   // alt: for (j=n-1;j>i;j--)
	aa[k][j] -= tmp*aa[i][j];
      }
      bb[k] -= tmp*bb[i];
    }
  }
  // Do back substitution
  for (i=n-1;i>=0;i--) {
    x[i]=bb[i];
    for (j=n-1;j>i;j--) {
      x[i] -= aa[i][j]*x[j];
    }
    x[i] /= aa[i][i];
  }
  // Deallocate memory
 _100:
  for (i=0;i<n;i++) {
    delete [] aa[i];
  }
  delete [] aa;
  delete [] bb;
  return retval;
}

int fsolve::broyden(void (*f)(double *x,double *fv,int n),
		    double *x0,double *f0,int n,double *eps,int *iter)
{
  double **A,*x1,*f1,*s,d,tmp;
  int i,j,k;
  int error_code = 1;
  
  // Allocate temporary memory
  x1 = new double [n];
  f1 = new double [n];
  s = new double [n];
  A = new double *[n];
  for (i=0;i<n;i++) {
    A[i] = new double [n];
  }
  
  // Create identity matrix for startup (consider Jacobian alternative)
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      A[i][j] = 0.0;
    }
    A[i][i] = 1.0;
  }
  // Main loop
  double n_sqr_eps = float(n) * *eps * *eps;
  for (k=0;k< *iter;k++){
    f(x0,f0,n);
    gelimd2(A,f0,s,n);      // Signs of 'f0', 's' are reversed
    d = 0.0;
    for (i=0;i<n;i++) {
      x1[i] = x0[i] - s[i];
      d += s[i]*s[i];
    }
    // If the RMS of the residuals is less than the solution tolerance, we're done
    if (d <= n_sqr_eps)
      {
	error_code = 0;
	break;
      }
    f(x1,f1,n);
    // Update A
    for (i=0;i<n;i++) {
      tmp = f1[i]+f0[i];
      for (j=0;j<n;j++) {
	tmp -= A[i][j]*s[j];
      }
      for (j=0;j<n;j++) {
	A[i][j] -= tmp*s[j]/d;
      }
    }
    for (i=0;i<n;i++) {
      f0[i] = f1[i];
      x0[i] = x1[i];
    }
  }
  // Replace maximum iterations with the number of iterations
  // required to find the solution
  *iter = k;
  
  // Deallocate memory
  for (i=0;i<n;i++) {
    delete [] A[i];
  }
  delete [] A;
  delete [] x1;
  delete [] f1;
  delete [] s;
  
  // Error codes are:
  // 0: successfully converged
  // 1: didn't converge within tolerance
  return error_code;
}
