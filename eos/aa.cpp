#include "aa.hpp"

/*
 *  Base class for the liquid metal EOS detailed in Anderson and Ahrens (1994).
 *  Slightly simplified in terms of C_v
 *  This equation of state is described by a fourth order BM isentropic EoS
 *  (V_0, KS, KS', KS''), a description of the volumetric heat capacity
 *  (which gives energy along an isochore as a function of temperature), and
 *  a description of the gruneisen parameter as a function of volume and energy
 *  (which defines the thermal pressure).
 */


double aa::gibbs(double pressure, double temperature)
{
  double volume = _volume(pressure, temperature);
  return helmholtz_free_energy(volume, temperature) + pressure * volume;
}

double aa::helmholtz_free_energy(double volume, double temperature)
{
  /*      
   * Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
   */
        
  return internal_energy(volume, temperature) - temperature*entropy(volume, temperature);
}

double aa::internal_energy(double volume, double temperature)
{
  // Returns the internal energy at the pressure and temperature of the mineral [J/mol]
  double Ts = _isentropic_temperature(volume);
  double E = (E_0 + _isentropic_energy_change(volume)
       + _isochoric_energy_change(Ts, temperature, volume));
            
  return E;
}

double aa::entropy(double volume, double temperature)
{
  // Returns the entropy at the pressure and temperature of the mineral [J/K/mol]
  double Ts = _isentropic_temperature(volume);
  double S = (S_0 + _entropy_kin(Ts, temperature, volume)
             + _entropy_el(Ts, temperature, volume)
	      + _entropy_pot(Ts, temperature, volume));

  return S;
}
      
// Potential heat capacity functions
aa::potential_properties aa::_lambdaxi(double V)
{
  potential_properties pot;
  
  double rhofrac = V_0/V;
        
  pot.xi = xi_0*std::pow(rhofrac, -0.6);
    //F = 1./(1. + np.exp((rhofrac - params['F'][0])/params['F'][1]))
    //lmda = (F*(params['lmda'][0]*rhofrac + params['lmda'][1]) + params['lmda'][2])*np.power(rhofrac, 0.4)
  pot.lmda=0.;
  return pot;
}

// Fourth order BM functions
aa::eos_properties aa::_rhofracxksis(double V)
{
  eos_properties eos;
  eos.rhofrac = V_0/V; // # rho/rho0 = V0/V
  eos.x = std::pow(eos.rhofrac, 1./3.); // equation 18
  eos.ksi1 = 0.75*(4. - Kprime_S); // equation 19
  eos.ksi2 = 0.375*(K_S*Kprime_prime_S
		       + Kprime_S*(Kprime_S - 7.)) + 143./24.; // # equation 20
  return eos;
}
         
// Contributions to the heat capacity

// High temperature limit of the kinetic contribution to the heat capacity
// Anderson and Ahrens (1994), just after equation 29.
double aa::_C_v_kin(double V, double T)
{
  return 1.5*n*constants::gas_constant;
}

// Equation A1
double aa::_C_v_el(double V, double T)
{
  return electronic::heat_capacity_v(T, V_0/V, T_el, Cv_el);
}

// Equation A15
double aa::_C_v_pot(double V, double T)
{
  potential_properties pot = _lambdaxi(V);
  double C_pot = (pot.lmda*T + pot.xi*theta) / (theta + T);
  return C_pot;
}

    
// Contributions to the internal energy
    
        
double aa::_internal_energy_kin(double Ts, double T, double V)
{
  double E_kin = 1.5*n*constants::gas_constant*(T - Ts);
  return E_kin;
}
    
double aa::_internal_energy_el(double Ts, double T, double V)
{
  return electronic::thermal_energy(T, V_0/V, T_el, Cv_el)	\
    - electronic::thermal_energy(Ts, V_0/V, T_el, Cv_el);
}
    
double aa::_internal_energy_pot(double Ts, double T, double V)
{
  potential_properties pot = _lambdaxi(V);
  double E_pot = (pot.lmda*(T - Ts) + theta*(pot.xi - pot.lmda)*std::log((theta + T)/(theta + Ts)));
  return E_pot;
}
    
// Contributions to entropy
double aa::_entropy_kin(double Ts, double T, double V)
{
  double S_kin = 0.;
  
  if (std::abs(T- Ts) > std::numeric_limits<double>::epsilon())
    {
      S_kin = 1.5*n*constants::gas_constant*(std::log(T) - std::log(Ts));
    }
  return S_kin;
}
        
double aa::_entropy_el(double Ts, double T, double V)
{    
  return electronic::entropy(T, V_0/V, T_el, Cv_el)	\
    - electronic::entropy(Ts, V_0/V, T_el, Cv_el);
}
    
double aa::_entropy_pot(double Ts, double T, double V)
{
  double S_pot = 0.;
    
  if (std::abs(T- Ts) > std::numeric_limits<double>::epsilon())
    {
      potential_properties pot = _lambdaxi(V);
      S_pot = (pot.lmda*std::log((theta + T)/(theta + Ts)) + pot.xi*std::log((T*(theta + Ts))/(Ts*(theta + T))));
    }
  return S_pot;
}
            
    
// Isentropic and isochoric calculations
    
// Temperature along an isentrope (Anderson and Ahrens; Equation B5)
double aa::_isentropic_temperature(double V)
{
  eos_properties eos = _rhofracxksis(V);
	  
  // equation B6 -- B10
  double a1 = eos.ksi2 / 8.;
  double a2 = ( eos.ksi1 + 3. * eos.ksi2 ) / 6.;
  double a3 = ( 1. + 2.*eos.ksi1 + 3.*eos.ksi2 ) / 4.;
  double a4 = (1. + eos.ksi1 + eos.ksi2)/2.;
  double a5 = (6. + 4.*eos.ksi1 + 3.*eos.ksi2)/24.;
    
  // equation B5
  double Ts = T_0*std::exp(grueneisen_0*std::log(eos.rhofrac)
			   + 13.5*grueneisen_prime*V_0*K_S *
			   (   (a1/(3*grueneisen_n + 8.))*(std::pow(eos.x,(3*grueneisen_n + 8.)) - 1.)
			       - (a2/(3*grueneisen_n + 6.))*(std::pow(eos.x,(3*grueneisen_n + 6.)) - 1.)
			       + (a3/(3*grueneisen_n + 4.))*(std::pow(eos.x,(3*grueneisen_n + 4.)) - 1.)
			       - (a4/(3*grueneisen_n + 2.))*(std::pow(eos.x,(3*grueneisen_n + 2.)) - 1.)
			       + (a5/(3*grueneisen_n + 0.))*(std::pow(eos.x,(3*grueneisen_n + 0.)) - 1.)));
                                
  return Ts;
}

// Pressure along the reference isentrope
double aa::_isentropic_pressure(double V)
{
  eos_properties eos = _rhofracxksis(V);
  double x2 = eos.x*eos.x;
  double x3 = eos.x*eos.x*eos.x;
  double x5 = x3*x2;
  double x7 = x5*x2;
  
  return 1.5*K_S * (x7 - x5) * (1. + eos.ksi1 - eos.ksi1*x2 + eos.ksi2 * (x2 - 1.) * (x2 - 1.));
}

// Birch Murnaghan equation of state expression for the energy change along an isentrope
// Anderson and Ahrens, 1994 (Equation 21)
double aa::_isentropic_energy_change(double V)
{
  eos_properties eos = _rhofracxksis(V);
  double x2 = eos.x*eos.x;
  double x4 = x2*x2;
  double x6 = x4*x2;
  double x8 = x4*x4;
        
  return 4.5*V_0*K_S * ((eos.ksi1 + 1.) * (x4/4. - x2/2. + 0.25) - eos.ksi1*(x6/6. - x4/4. + 1./12.)
			+ eos.ksi2*(x8/8. - x6/2. + 0.75*x4 - x2/2. + 0.125));
}

// int Cv dT 
double aa::_isochoric_energy_change(double Ts, double T, double V)
{
  return _internal_energy_kin(Ts, T, V)
    + _internal_energy_el(Ts, T, V)
    + _internal_energy_pot(Ts, T, V);
}
    


double aa::_volume(double pressure, double temperature)
{
  // Returns molar volume. :math:`[m^3]`
  return _find_volume(0.1*V_0, 10.*V_0, 1.e-12, pressure, temperature);
}

double aa::_pressure(double volume, double temperature)
{
  // Returns the pressure of the mineral at a given temperature and volume [Pa]
       
  double Ts = _isentropic_temperature(volume);
  double dE = _isochoric_energy_change(Ts, temperature, volume);
  double E1 = _isentropic_energy_change(volume); // should also include E_0 given the expression in Anderson and Ahrens. Here, we take the energy change relative to the reference isentrope (effective E_0 = 0). The energy at standard state is *only* used to calculate the final energies, not the physical properties.
  double E2 = E1 + dE;
  double dP = (grueneisen_0*dE + 0.5*grueneisen_prime*std::pow(V_0/volume,grueneisen_n)*(E2*E2 - E1*E1))/volume;
  return _isentropic_pressure(volume) + dP;
}


double aa::_find_volume(double a, double b, double t, double pressure, double temperature)
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






	      



  
