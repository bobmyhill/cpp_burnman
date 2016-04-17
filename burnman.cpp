#include "burnman.hpp"

int main ()
{

  slb hcp_iron = models::hcp_iron_model();
  aa liq_iron = models::liq_iron_model();


  std::cout << hcp_iron.name << std::endl;
  std::cout << hcp_iron.gibbs(1.e5, 1809.) << std::endl;
  std::cout << hcp_iron.gibbs(1.e9, 2000.) << std::endl;
  std::cout << hcp_iron.gibbs(1.e9, 3000.) << std::endl;
  std::cout << hcp_iron.gibbs(1.e10, 3000.) << std::endl;
  std::cout << liq_iron.name << std::endl;
  std::cout << liq_iron.gibbs(1.e5, 1809.) << std::endl;
  std::cout << liq_iron.gibbs(1.e9, 2000.) << std::endl;
  std::cout << liq_iron.gibbs(1.e9, 3000.) << std::endl;
  std::cout << liq_iron.gibbs(1.e10, 3000.) << std::endl;

  
  extended_margules fesio = models::fesio_model();
  boost::numeric::ublas::vector<double> molar_fractions(3);
  molar_fractions(0) = 0.1;
  molar_fractions(1) = 0.7;
  molar_fractions(2) = 0.2;
  
    
  boost::numeric::ublas::vector<double> a = fesio.activities(1.e10, 3000., molar_fractions);

  

  return 0;
}


