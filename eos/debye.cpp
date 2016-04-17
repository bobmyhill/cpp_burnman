#include "debye.hpp"
// Functions for the Debye model of a solid.
#include <iostream>

namespace debye {
  
  double _chebval(double x, double *c, int csize);
  double debye_fn_cheb(double x);
  double thermal_energy(double T, double debye_T, double n);
  double entropy(double T, double debye_T, double n);
  double helmholtz_free_energy(double T, double debye_T, double n);
  double heat_capacity_v(double T, double debye_T, double n);


  double chebyshev_representation[17] = {2.707737068327440945/2.0, 0.340068135211091751, -0.12945150184440869e-01,
					 0.7963755380173816e-03, -0.546360009590824e-04, 0.39243019598805e-05, 
					 -0.2894032823539e-06, 0.217317613962e-07, -0.16542099950e-08, 
					 0.1272796189e-09, -0.987963460e-11, 0.7725074e-12, -0.607797e-13, 
					 0.48076e-14, -0.3820e-15, 0.305e-16, -0.24e-17};
};


double debye::_chebval(double x, double *c, int csize)
{
  /*
   * Evaluate a Chebyshev series at points x.
   */

  double c0, c1, tmp, x2;

  if (csize == 1)
    {
      c0 = c[0];
      c1 = 0;
    }
  else if (csize == 2)
    {
      c0 = c[0];
      c1 = c[1];
    }
  else
    {
      x2 = 2*x;
      c0 = c[csize-2];
      c1 = c[csize-1];
      for (int i=3; i < csize + 1; i++)
	{
	  tmp = c0;
	  c0 = c[csize-i] - c1;
	  c1 = tmp + c1*x2;
	}
    }
  return c0 + c1*x;
}


double debye::debye_fn_cheb(double x)
{
  /*
   * Evaluate the Debye function using a Chebyshev series expansion coupled with
   * asymptotic solutions of the function.  Shamelessly adapted from the GSL implementation
   * of the same function (Itself adapted from Collected Algorithms from ACM).
   * Should give the same result as debye_fn(x) to near machine-precision.
   */

  double eps = std::numeric_limits<double>::epsilon();
  double sqrt_eps = std::sqrt(eps);
  double log_eps = std::log(eps);
  
  double val_infinity = 19.4818182068004875;
  double xcut = -log_eps;

  double f;
  
  assert(x > 0.0); // check for invalid x

  if (x < 2.0 * std::sqrt(2.0) * sqrt_eps)
    {
      f = 1.0 - 3.0*x/8.0 + x*x/20.0;
    }
  else if (x <= 4.0)
    {
      double t = x*x/8.0 - 1.0;
      int csize = sizeof(chebyshev_representation)/sizeof(*chebyshev_representation);
      double c = _chebval(t, chebyshev_representation, csize);
      f = c - 0.375*x;
    }
  else if (x < -(std::log(2.0) + log_eps ))
    {
      int nexp = int(std::floor(xcut/x));
      double ex  = std::exp(-x);
      double xk  = nexp * x;
      double rk  = nexp;
      double sum = 0.0;

      double xk_inv;
      for (int i=nexp; i>0; i-=1)
	{
	  xk_inv = 1.0/xk;
	  sum *= ex;
	  sum += (((6.0*xk_inv + 6.0)*xk_inv + 3.0)*xk_inv + 1.0) / rk;
	  rk -= 1.0;
	  xk -= x;
	}
      f = val_infinity/(x*x*x) - 3.0 * sum * ex;
    }	
  else if (x < xcut)
    {
      double x3 = x*x*x;
      double sum = 6.0 + 6.0*x + 3.0*x*x + x3;
      f = (val_infinity - 3.0 * sum * std::exp(-x)) / x3;
    }
  else
    {
      f = ((val_infinity/x)/x)/x;
    }

  return f;
}


double debye::thermal_energy(double T, double debye_T, double n)
{
  /*
   * calculate the thermal energy of a substance.  Takes the temperature,
   * the Debye temperature, and n, the number of atoms per molecule.
   * Returns thermal energy in J/mol
   */
  double E_th = 0.;
  if (T > std::numeric_limits<double>::epsilon())
    {
      E_th = 3.*n*constants::gas_constant*T * debye_fn_cheb(debye_T/T);
    }
  return E_th;
}

double debye::heat_capacity_v(double T,double debye_T,double n)
{
  /* 
   * Heat capacity at constant volume.  In J/K/mol
   */
  double C_v = 0.;
  if (T > std::numeric_limits<double>::epsilon())
    {
      double x = debye_T/T;
      C_v = 3.0*n*constants::gas_constant* ( 4.0*debye_fn_cheb(x) - 3.0*x/(std::exp(x)-1.0) );
    }
  
  return C_v;
}

double debye::helmholtz_free_energy(double T, double debye_T, double n)
{
  /*
   * Helmholtz free energy of lattice vibrations in the Debye model.
   * It is important to note that this does NOT include the zero 
   * point energy of vibration for the lattice.  As long as you are 
   * calculating relative differences in F, this should cancel anyways.
   * In Joules.
   */
  
  double F = 0.;
      
  if (T > std::numeric_limits<double>::epsilon())
    {
      double x = debye_T/T;
      F = n * constants::gas_constant * T * ( 3.0 * std::log( 1.0 - std::exp(-x)) - debye_fn_cheb(x) );
    }
  return F;  
}

double debye::entropy(double T, double debye_T, double n)
{
  /*
   * Entropy due to lattice vibrations in the Debye model [J/K]
   */
  
  double S = 0.;
  
  if (T > std::numeric_limits<double>::epsilon())
    {
      double x = debye_T/T;
      S = n * constants::gas_constant * ( 4. * debye_fn_cheb(x) - 3. * std::log( 1.0 - std::exp(-x) ) );
    }
  
  return S;
}

