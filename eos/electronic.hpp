#include <cmath>

namespace electronic {
  double thermal_energy (double T, double x, double Tel_0, double Cvel_max);
  double entropy(double T, double x, double Tel_0, double Cvel_max);
  double helmholtz_free_energy(double T, double x, double Tel_0, double Cvel_max);
  double heat_capacity_v(double T, double x, double Tel_0, double Cvel_max);

};
