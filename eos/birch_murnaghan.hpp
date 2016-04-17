#ifndef __burnman__bm_h
#define __burnman__bm_h

#include <cmath>
#include <sstream>
#include <algorithm>

class bm3 {
public:
  std::string name;
  double P_0;
  double V_0;
  double K_0;
  double Kprime_0;
  
  double isothermal_bulk_modulus(double volume);
  double birch_murnaghan(double x);
  double volume(double pressure);
  double pressure(double volume);
  double _find_x(double a, double b, double t, double pressure);

};

#endif
