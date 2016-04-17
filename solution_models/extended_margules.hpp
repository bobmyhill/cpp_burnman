#ifndef __burnman__extended_margules_h
#define __burnman__extended_margules_h

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "../global.hpp"
#include <string>


class extended_margules {
public:
  std::string name;

  int n_endmembers;
  int n_occupancies;

  boost::numeric::ublas::matrix<double> endmember_occupancies; // size n_endmembers, n_occupancies
  boost::numeric::ublas::vector<double> site_multiplicities; // size n_occupancies
  
  boost::numeric::ublas::matrix<double> We; // size n_endmembers, n_endmembers
  boost::numeric::ublas::matrix<double> Ws; // size n_endmembers, n_endmembers
  boost::numeric::ublas::matrix<double> Wv; // size n_endmembers, n_endmembers
  boost::numeric::ublas::matrix<double> Wk; // size n_endmembers, n_endmembers


  boost::numeric::ublas::vector<double> _non_ideal_function(boost::numeric::ublas::matrix <double> W, boost::numeric::ublas::vector<double> molar_fractions);
  boost::numeric::ublas::vector<double> activity_coefficients(double pressure, double temperature, boost::numeric::ublas::vector<double> molar_fractions);
  boost::numeric::ublas::vector<double> activities(double pressure, double temperature, boost::numeric::ublas::vector<double> molar_fractions);
  boost::numeric::ublas::vector<double> ideal_activities(boost::numeric::ublas::vector<double> molar_fractions);
  boost::numeric::ublas::vector<double> _calculate_endmember_configurational_entropies();
};



extended_margules extended_margules_model(std::string name, int n_endmembers, int n_occupancies,
					  double *occ, double *mult,
					  double *energy, double *entropy, double *volume, double *modulus);

#endif
