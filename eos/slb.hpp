#ifndef __burnman__slb_h
#define __burnman__slb_h

#include <string>
using namespace std;

class slb {
public:
  string name;
  double P_0;
  double T_0;
  double F_0;
  double V_0;
  double K_0;
  double Kprime_0;
  double Debye_0;
  double grueneisen_0;
  double q_0;
  double Cv_el;
  double T_el;
  double n;
  double molar_mass;
  
  double gibbs (double pressure, double temperature);
  double helmholtz_free_energy(double volume, double temperature);
  double _pressure(double volume, double temperature);
  double _volume(double pressure, double temperature);
  double _find_volume(double a, double b, double t, double pressure, double temperature);
  double _debye_temperature(double x);
};

#endif
