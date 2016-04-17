#include "slb.hpp"
#include "electronic.hpp"
#include "debye.hpp"
#include "birch_murnaghan.hpp"
#include <iostream>

double slb::gibbs(double pressure, double temperature)
{
  double volume = _volume(pressure, temperature);
  double F = helmholtz_free_energy(volume, temperature);
  return F + pressure * volume;
}


double slb::helmholtz_free_energy(double volume, double temperature)
{
  // Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
  
  double x = V_0 / volume;
  double f = 1./2. * (std::pow(x, 2./3.) - 1.);
  double Debye_T = _debye_temperature(V_0/volume);

  double F_quasiharmonic = (debye::helmholtz_free_energy(temperature, Debye_T, n) 
			    + electronic::helmholtz_free_energy(temperature, V_0/volume,
								T_el, Cv_el)) 
    - (debye::helmholtz_free_energy(T_0, Debye_T, n)		
       + electronic::helmholtz_free_energy(T_0, V_0/volume,
					   T_el, Cv_el));
        
  double b_iikk= 9.*K_0;// # EQ 28
  double b_iikkmm= 27.*K_0*(Kprime_0-4.);// # EQ 29

  double F = F_0 +							\
    0.5*b_iikk*f*f*V_0 + (1./6.)*V_0*b_iikkmm*f*f*f +			\
    F_quasiharmonic;

  return F;
}

double slb::_pressure(double volume, double temperature)
{
  // Numerical differentiation. Should be reasonably accurate.
  double dV = V_0*1.e-3;
  double F0 = helmholtz_free_energy(volume-0.5*dV, temperature);
  double F1 = helmholtz_free_energy(volume+0.5*dV, temperature);
  return (F0 - F1)/dV;
}


double slb::_volume(double pressure, double temperature)
{
  return _find_volume(0.1*V_0, 10.*V_0, 1.e-12, pressure, temperature);
}

double slb::_find_volume(double a, double b, double t, double pressure, double temperature)
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
  fa = _pressure(sa, temperature) - pressure;
  fb = _pressure(sb, temperature) - pressure;
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

    fb = _pressure(sb, temperature) - pressure;
    
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

double slb::_debye_temperature(double x)
{
  /*
   *    Finite strain approximation for Debye Temperature [K]
   *    x = ref_vol/vol
   */
  double f = 1./2. * (std::pow(x, 2./3.) - 1.);
  double a1_ii = 6. * grueneisen_0; //  EQ 47
  double a2_iikk = -12.*grueneisen_0+36.*std::pow(grueneisen_0,2.) - 18.*q_0*grueneisen_0; //EQ 47
  return Debye_0 * std::sqrt(1. + a1_ii * f + 1./2. * a2_iikk*f*f);
}
