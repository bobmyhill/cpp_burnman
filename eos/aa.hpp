#ifndef __burnman__aa_h
#define __burnman__aa_h

#include <string>
#include <cmath>
#include <limits>
#include "electronic.hpp"
#include "../global.hpp"

class aa {
public:
  std::string name;
  double P_0;
  double T_0;
  double S_0;
  double molar_mass;
  
  double V_0;
  double E_0;
  double K_S;
  double Kprime_S;
  double Kprime_prime_S;
  double grueneisen_0;
  double grueneisen_prime;
  double grueneisen_n;
  double T_el;
  double Cv_el;
  double theta;
  double xi_0;
  double n;

  struct eos_properties {
    double rhofrac;
    double x;
    double ksi1;
    double ksi2;
  };
  
  struct potential_properties {
    double lmda;
    double xi;
  };
  
  double gibbs (double pressure, double temperature);
  double helmholtz_free_energy(double volume, double temperature);
  double internal_energy(double volume, double temperature);
  double entropy(double volume, double temperature);
  potential_properties _lambdaxi(double V);
  eos_properties _rhofracxksis(double V);
  double _C_v_kin(double V, double T);
  double _C_v_el(double V, double T); 
  double _C_v_pot(double V, double T);
  double _internal_energy_kin(double Ts, double T, double V);
  double _internal_energy_el(double Ts, double T, double V);
  double _internal_energy_pot(double Ts, double T, double V);
  double _entropy_kin(double Ts, double T, double V);
  double _entropy_el(double Ts, double T, double V);
  double _entropy_pot(double Ts, double T, double V);
  double _isentropic_temperature(double V);
  double _isentropic_pressure(double V);
  double _isentropic_energy_change(double V);
  double _isochoric_energy_change(double Ts, double T, double V);
  double _volume(double pressure, double temperature);
  double _pressure(double volume, double temperature);
  double _find_volume(double a, double b, double t, double pressure, double temperature);
};

#endif
