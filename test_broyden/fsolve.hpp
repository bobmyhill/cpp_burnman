#include <cmath>

namespace fsolve {
  int gelimd2(double **a,double *b,double *x, int n);
  
  int broyden(void (*f)(double *x,double *fv,int n),double *x,
	      double *fv,int n,double *eps,int *iter);
}
