#include <cmath>
#include "../global.hpp"

namespace einstein {
  double thermal_energy(double T, double einstein_T, double n);
  double entropy(double T, double einstein_T, double n);
  double heat_capacity_v(double T, double einstein_T, double n);
};
