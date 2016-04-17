#include "einstein.hpp"

// Functions for the Einstein model of a solid.

namespace einstein {
  double thermal_energy(double T, double einstein_T, double n);
  double entropy(double T, double einstein_T, double n);
  double heat_capacity_v(double T, double einstein_T, double n);
  
  double eps = 1.e-12;
};

double einstein::thermal_energy(double T, double einstein_T, double n)
{
    
  // calculate the thermal energy of a substance.  Takes the temperature,
  // the Einstein temperature, and n, the number of atoms per molecule.
  // Returns thermal energy in J/mol

  double E_th = 3.*n*constants::gas_constant*einstein_T*0.5; // zero point energy
  
  if (T > eps)
    {
      double x = einstein_T/T;
      E_th = 3.*n*constants::gas_constant*einstein_T*( 0.5 + 1. / (std::exp( x ) - 1.0) ); // include the zero point energy
    }
  return E_th;
}

double einstein::entropy(double T, double einstein_T, double n)
{
  //   Entropy at constant volume.  In J/K/mol

  double S = 0.;
  if (T > eps)
    {
      double x = einstein_T/T;
      S = 3.0*n*constants::gas_constant * (x*(1./(std::exp(x) - 1) + 1.) - std::log(std::exp(x) - 1.));
    }
  return S;
}

double einstein::heat_capacity_v(double T,double einstein_T, double n)
{
  //    Heat capacity at constant volume.  In J/K/mol
  double C_v = 0.; 
  if (T > eps)
    {
      double x = einstein_T/T;
      C_v = 3.0*n*constants::gas_constant* ( x * x * std::exp( x ) / std::pow( std::exp( x ) - 1.0, 2.0 ) );
    }
  return C_v;
}

