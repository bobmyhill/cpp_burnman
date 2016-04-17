#include "extended_margules.hpp"
#include <iostream>
#include <sstream>

boost::numeric::ublas::vector<double> extended_margules::_non_ideal_function(boost::numeric::ublas::matrix <double> W,
									     boost::numeric::ublas::vector<double> molar_fractions)
{
  // equation (6') of Helffrich and Wood, 1989
  boost::numeric::ublas::zero_vector<double> zero (n_endmembers);
  boost::numeric::ublas::vector<double> RTlny = zero;
  double val;

  for (int l=0; l<n_endmembers; l++)
    {
      val = 0.;
      for (int i=0; i<n_endmembers; i++)
	{
	  if (i != l)
	    {
	      val += 0.5 * molar_fractions(i) * (W(l,i) * (1. - molar_fractions(l) + molar_fractions(i) + 2. * molar_fractions(l) * (molar_fractions(l) - molar_fractions(i) - 1.)) + W(i,l) * (1. - molar_fractions(l) - molar_fractions(i) - 2. * molar_fractions(l) * (molar_fractions(l) - molar_fractions(i) - 1.)));
	    }
	  for (int j=i+1; j<n_endmembers; j++)
	    {
	      if (j != l)
		{
		  val += molar_fractions(i) * molar_fractions(j) * (W(i,j) * (molar_fractions(i) - molar_fractions(j) - 0.5) + W(j,i) * (molar_fractions(j) - molar_fractions(i) - 0.5));
		}
	    }
	}
      RTlny(l) = val;
    }
  return RTlny;
}

  
boost::numeric::ublas::vector<double> extended_margules::activity_coefficients(double pressure, double temperature,
									       boost::numeric::ublas::vector<double> molar_fractions)
{
  
  double Kprime = 4.;

  boost::numeric::ublas::matrix <double> intVdP (n_endmembers, n_endmembers);
  boost::numeric::ublas::matrix <double> Wg (n_endmembers, n_endmembers);

  for (int i=0; i<n_endmembers; i++)
    {
      for (int j=0; j<n_endmembers; j++)
	{
	  intVdP(i,j) = Wv(i,j)*Wk(i,j)*(1. - std::pow(1. + Kprime/Wk(i,j)*pressure, (1. - 1./Kprime)))/(1. - Kprime);
	  Wg(i,j) = We(i,j) - temperature*Ws(i,j) + intVdP(i,j);
	}
    }
  
  boost::numeric::ublas::vector<double> partial_nonideal_excess_gibbs = _non_ideal_function(Wg, molar_fractions);

  
  boost::numeric::ublas::vector<double> gammas (n_endmembers);
  for (int i=0; i<n_endmembers; i++)
    {
      gammas(i) =  std::exp(partial_nonideal_excess_gibbs(i) /
			    8.31446*temperature);
    }
  
  return gammas;
}

boost::numeric::ublas::vector<double> extended_margules::ideal_activities(boost::numeric::ublas::vector<double> molar_fractions)
{
  boost::numeric::ublas::vector<double> endmember_configurational_entropies = _calculate_endmember_configurational_entropies();
  boost::numeric::ublas::vector<double> site_occupancies = prod(molar_fractions, endmember_occupancies);
  boost::numeric::ublas::vector<double> ideal_activities (n_endmembers);

  for (int e = 0; e < n_endmembers; e++)
    {
      ideal_activities[e] = 1.0;
      for (int occ = 0; occ < n_occupancies; occ++)
	{
	  if (endmember_occupancies(e,occ) > 1e-10)
	    {
	      ideal_activities[e] = ideal_activities(e)			\
		* std::pow(site_occupancies(occ), endmember_occupancies(e,occ) \
			   * site_multiplicities(occ));
	    }
	}
      ideal_activities[e] = std::exp(endmember_configurational_entropies(e) / 8.31446) * ideal_activities(e);
    }
  return ideal_activities;

}

boost::numeric::ublas::vector<double> extended_margules::activities(double pressure, double temperature, boost::numeric::ublas::vector<double> molar_fractions)
{
  
  boost::numeric::ublas::vector<double> a (n_endmembers);
  boost::numeric::ublas::vector<double> ideal = ideal_activities(molar_fractions);
  boost::numeric::ublas::vector<double> gamma = activity_coefficients(pressure, temperature, molar_fractions);
  for  (int i = 0; i < n_endmembers; i++)
    {
      a(i) = ideal(i)*gamma(i);
    }
  
  return a;
}

boost::numeric::ublas::vector<double> extended_margules::_calculate_endmember_configurational_entropies()
{
  boost::numeric::ublas::vector<double> endmember_configurational_entropies (n_endmembers);
    
  for (int i = 0; i < n_endmembers; i++)
    {
      endmember_configurational_entropies(i) = 0;
      for (int occ = 0; occ < n_occupancies; occ++)
	{
	  if (endmember_occupancies(i,occ) > 1e-10)
	    {
	      endmember_configurational_entropies(i) -= 		\
		8.31446 * site_multiplicities(occ)	\
		* endmember_occupancies(i,occ) * std::log(endmember_occupancies(i,occ));
	    }
	}
    }

  return endmember_configurational_entropies;
}



extended_margules extended_margules_model(std::string name, int n_endmembers, int n_occupancies,
					  double *occ, double *mult,
					  double *energy, double *entropy, double *volume, double *modulus)
{
  extended_margules model;
  model.name=name;
  model.n_endmembers = n_endmembers;
  model.n_occupancies = n_occupancies;
  
  model.endmember_occupancies.resize(n_endmembers, n_occupancies);
  for (int i=0; i<n_endmembers; i++)
    {
      for (int j=0; j<n_occupancies; j++)
	{
	  model.endmember_occupancies(i,j) = occ[i*n_occupancies + j]; 
	}
    }
    
  model.site_multiplicities.resize(n_occupancies);
  for (int i=0; i<n_occupancies; i++)
    {
      model.site_multiplicities(i) = mult[i];
    }

  model.We.resize(n_endmembers, n_endmembers);
  model.Ws.resize(n_endmembers, n_endmembers);
  model.Wv.resize(n_endmembers, n_endmembers);
  model.Wk.resize(n_endmembers, n_endmembers);
  
  for (int i=0; i<n_endmembers; i++)
    {
      for (int j=0; j<n_endmembers; j++)
	{
	  model.We(i,j) = energy[i*n_endmembers + j]; 
	  model.Ws(i,j) = entropy[i*n_endmembers + j]; 
	  model.Wv(i,j) = volume[i*n_endmembers + j]; 
	  model.Wk(i,j) = modulus[i*n_endmembers + j]; 
	}
    }
  
  return model;
}

