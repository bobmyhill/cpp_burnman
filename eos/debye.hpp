#include <cmath>
#include <assert.h>
#include <limits>
#include <vector>
#include "../global.hpp"

namespace debye {
  double _chebval(double x, double *c, int csize);
  double debye_fn_cheb(double x);
  double thermal_energy(double T, double debye_T, double n);
  double entropy(double T, double debye_T, double n);
  double helmholtz_free_energy(double T, double debye_T, double n);
  double heat_capacity_v(double T, double debye_T, double n);
  
};
