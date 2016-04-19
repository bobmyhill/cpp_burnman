#include <cmath>
#include <iostream>
#include <limits>

#include "fsolve.hpp"

// x is the solution vector
// fv is the function vector
// n is the number of equations
void func(double *x,double *fv,int n)
{
  fv[0] = cos(2.0*x[0]) - cos(2.0*x[1]) -0.4;
  fv[1] = 2.0*(x[1]-x[0]) + sin(2.0*x[1]) - sin(2.0*x[0]) - 1.2;
}



int main()
{
  int iterations = 50; // maximum number of iterations
  int error_code;
  
  int n_equations = 2;
  
  double tol = 1e-10;
  double guesses[n_equations];
  double residuals[n_equations];
  
  guesses[0] = 5.;
  guesses[1] = 5.;
  error_code = fsolve::broyden(func,guesses,residuals,n_equations,&tol,&iterations);
  std::cout << "Tolerance: " << tol << std::endl;
  std::cout << "Error code: " << error_code << std::endl;
  std::cout << "Iterations: " << iterations << std::endl;
  std::cout << "x[0] = " << guesses[0] << std::endl;
  std::cout << "x[1] = " << guesses[1] << std::endl;
  std::cout << "r[0] = " << residuals[0] << std::endl;
  std::cout << "r[1] = " << residuals[1] << std::endl;
  return 0;
}
