#include "birch_murnaghan.hpp"
double bm3::isothermal_bulk_modulus(double volume)
{

  /* 
   * compute the bulk modulus as per the third order
   * birch-murnaghan equation of state.  Returns bulk
   * modulus in the same units as the reference bulk
   * modulus.  Pressure must be in :math:`[Pa]`.
   */
  
  double x = V_0/volume;
  double f = 0.5*(std::pow(x, 2./3.) - 1.0);
  
  double K = std::pow(1. + 2.*f, 5./2.) * (K_0 + (3. * K_0 * Kprime_0 - 5.*K_0 ) * f + 27./2. * (K_0*Kprime_0 - 4.*K_0)*f*f);
  return K;
  
}

double bm3::birch_murnaghan(double x)
{
  /*
   * equation for the third order birch-murnaghan equation of state, returns
   * pressure in the same units that are supplied for the reference bulk
   * modulus (params['K_0'])
   */
  
  return 3.*K_0/2. * (std::pow(x, 7./3.) - std::pow(x, 5./3.))		\
    * (1. - .75*(4.-Kprime_0 )*(std::pow(x, 2./3.) - 1.)) + P_0;
}

double bm3::volume(double pressure)
{
  /*
   * Get the birch-murnaghan volume at a reference temperature for a given
   * pressure :math:`[Pa]`. Returns molar volume in :math:`[m^3]`
   */

  //sol = bracket(delta_pressure, V_0, 1.e-2*V_0);
  //return opt.brentq(func, sol[0], sol[1])
  
  return V_0/_find_x(0.1, 10., 1.e-12, pressure);
}

double bm3::_find_x(double a, double b, double t, double pressure)
{
  double c;
  double d;
  double e;
  double fa;
  double fb;
  double fc;
  double m;
  double macheps = std::numeric_limits<double>::epsilon();
  double p;
  double q;
  double r;
  double s;
  double sa;
  double sb;
  double tol;
//
//  Make local copies of A and B.
//
  sa = a;
  sb = b;
  fa = birch_murnaghan(sa) - pressure;
  fb = birch_murnaghan(sb) - pressure;

  c = sa;
  fc = fa;
  e = sb - sa;
  d = e;


  for ( ; ; )
  {
    if ( std::abs ( fc ) < std::abs ( fb ) )
    {
      sa = sb;
      sb = c;
      c = sa;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    tol = 2.0 * macheps * std::abs ( sb ) + t;
    m = 0.5 * ( c - sb );

    if ( std::abs ( m ) <= tol || fb == 0.0 )
    {
      break;
    }

    if ( std::abs ( e ) < tol || std::abs ( fa ) <= std::abs ( fb ) )
    {
      e = m;
      d = e;
    }
    else
    {
      s = fb / fa;

      if ( sa == c )
      {
        p = 2.0 * m * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
        q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
      }

      if ( 0.0 < p )
      {
        q = - q;
      }
      else
      {
        p = - p;
      }

      s = e;
      e = d;

      if ( 2.0 * p < 3.0 * m * q - std::abs ( tol * q ) &&
        p < std::abs ( 0.5 * s * q ) )
      {
        d = p / q;
      }
      else
      {
        e = m;
        d = e;
      }
    }
    sa = sb;
    fa = fb;

    if ( tol < std::abs ( d ) )
    {
      sb = sb + d;
    }
    else if ( 0.0 < m )
    {
      sb = sb + tol;
    }
    else
    {
      sb = sb - tol;
    }

    fb = birch_murnaghan(sb) - pressure;
    if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
    {
      c = sa;
      fc = fa;
      e = sb - sa;
      d = e;
    }
  }
  return sb;
}
