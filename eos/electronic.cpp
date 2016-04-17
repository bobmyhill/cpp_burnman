#include "einstein.hpp"

// Functions for the electronic contribution to the free energy.
// The function is based on a semiempirical representation of the 
// electronic contribution to the electronic heat capacity. 

// The heat_capacity_v curve is an einstein function at high temperature, 
// reaching a user-defined maximum (as a multiple of gas_constant)
// at very high temperatures. This value should be ~2.5 - 3. 
// The temperature scaling is governed by a parameter Tel, where 
// Tel = Tel_0*(V_0/V)^(3/2). Larger compressions result in larger Tel 
// and thus smaller electronic contributions to the heat capacity, 
// in agreement with observations.


// At ~0.38*Tel, the tangent to the heat capacity curve passes 
// through the origin. At lower temperatures, the heat 
// capacity is approximated by this tangent. This satisfies the 
// requirement of XXXX that Cv_el be linear at low temperatures.


namespace electronic {
  double thermal_energy (double T, double x, double Tel_0, double Cvel_max);
  double entropy(double T, double x, double Tel_0, double Cvel_max);
  double helmholtz_free_energy(double T, double x, double Tel_0, double Cvel_max);
  double heat_capacity_v(double T, double x, double Tel_0, double Cvel_max);

  // the fraction of Tel below which the Cv curve is linear
  // this value originates from the form of the Einstein heat capacity.
  double f_linear = 0.388247;
  double eps = 1.e-12;
  
};



double electronic::thermal_energy (double T, double x, double Tel_0, double Cvel_max)
{
  // Calculate the thermal energy of a substance.  Takes the temperature,
  // the characteristic electronic temperature and the maximum heat capacity.
  // Returns thermal energy in J/mol

  // x = ref_vol/vol
  float E_th;
  float Tel = Tel_0 * std::pow(x, 3./2.);
  if (T < f_linear*Tel)
    {
      E_th = 0.5 * (38./3.*Cvel_max)/Tel * std::pow(T, 2.);
    }
  else
    {
      E_th = einstein::thermal_energy(T, Tel, Cvel_max/3.)		\
	- einstein::thermal_energy(f_linear*Tel, Tel, Cvel_max/3.)	\
	+ 0.5 * (38./3.*Cvel_max)/Tel * std::pow(f_linear*Tel, 2.) ;
    }
  // note that subtracting the low temperature part of the curve
  // also removes the unwanted zero point energy.
        
  return E_th;
}
        
double electronic::entropy(double T, double x, double Tel_0, double Cvel_max)
{
  // Entropy due to electronic density of states [J/K/mol]

  // x = ref_vol/vol
  float S;
  float Tel = Tel_0*std::pow(x, 3./2.);
  if (T < f_linear*Tel)
    {
      S = (38./3.*Cvel_max)/Tel * T;
    }
  else
    {
      S = einstein::entropy(T, Tel, Cvel_max/3.)	       \
	- einstein::entropy(f_linear*Tel, Tel, Cvel_max/3.)    \
	+ (38./3.*Cvel_max)/Tel * f_linear*Tel;
    }
  return S;
}

double electronic::helmholtz_free_energy(double T, double x, double Tel_0, double Cvel_max)
{
  // Helmholtz free energy from the electronic density of states.
  // It is important to note that this does NOT include the zero 
  // point energy.  As long as you are calculating relative differences 
  // in F, this should cancel anyways.
  // In Joules.

  float F;
  if (T <= eps)
    {
      F = 0.;
    }
  else
    {
      F = thermal_energy(T, x, Tel_0, Cvel_max) \
        - T*entropy(T, x, Tel_0, Cvel_max);
    }
    
  return F;
}

double electronic::heat_capacity_v(double T, double x, double Tel_0, double Cvel_max)
{
  // Heat capacity contribution from the electronic 
  // density of states at constant volume.  In J/K/mol

  // x = ref_vol/vol
    
  double Tel = Tel_0*std::pow(x, 3./2.);
  double C_v;
  if (T < f_linear*Tel)
    {
      C_v = (38./3.*Cvel_max)/Tel * T; // 38./3. includes gas_constant
    }
  else
    {
      C_v = einstein::heat_capacity_v(T, Tel, Cvel_max/3.);
    }

  return C_v;
}
