#include "models.hpp"

namespace models
{
  extended_margules fesio_model()
  {
    std::string name = "Fe-Si-O model";
    int n_endmembers = 3;
    int n_occupancies = 6;
    double occ [18] = {1., 0., 0., 1., 0., 0., // Fe0.5Fe0.5
		       1., 0., 0., 0., 1., 0., //Fe0.5Si0.5 (ordered)
		       0., 0., 1., 0., 0., 1.}; //FeO0.5FeO0.5
    double mult [6] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    
    double energy [9] = {0., 0., 0.,
			 0., 0., 0.,
			 0., 0., 0.};
    double entropy [9] = {0., 0., 0.,
			  0., 0., 0.,
			  0., 0., 0.};
    double volume [9] = {0., 0., 0.,
			 0., 0., 0.,
			 0., 0., 0.};
    double modulus [9] = {100.e9, 100.e9, 100.e9,
			  100.e9, 100.e9, 100.e9,
			  100.e9, 100.e9, 100.e9};
    
    extended_margules fesio = extended_margules_model(name, n_endmembers, n_occupancies, occ, mult, energy, entropy, volume, modulus);
    return fesio;
  }
  
  
  

  slb hcp_iron_model()
  {
    slb hcp_iron;
    hcp_iron.name = "hcp iron";
    hcp_iron.P_0 = 1.e5;
    hcp_iron.T_0 = 300.;
    hcp_iron.F_0 = -4165.;
    hcp_iron.V_0 = 6.764e-6;
    hcp_iron.K_0 = 161.9e9;
    hcp_iron.Kprime_0 = 5.15;
    hcp_iron.Debye_0 = 395.;
    hcp_iron.grueneisen_0 = 2.0;
    hcp_iron.q_0 = 1.0;
    hcp_iron.Cv_el = 3.0;
    hcp_iron.T_el = 6000.;
    hcp_iron.n = 1.;
    hcp_iron.molar_mass = 0.055845;
    
    return hcp_iron;
  }
  
  aa liq_iron_model()
  {
    aa liq_iron;
    
    liq_iron.name = "liquid iron";
    liq_iron.P_0 = 1.e5;
    liq_iron.T_0 = 1809.;
    liq_iron.S_0 = 99.823; 
    liq_iron.molar_mass = 0.055845;
    liq_iron.V_0 = 0.055845/7019.;
    liq_iron.E_0 = 72700.;
    liq_iron.K_S = 109.7e9;
    liq_iron.Kprime_S = 4.661;
    liq_iron.Kprime_prime_S = -4.661/109.7e9;
    liq_iron.grueneisen_0 = 1.735;
    liq_iron.grueneisen_prime = -0.130/0.055845*1.e-6;
    liq_iron.grueneisen_n = -1.870;
    liq_iron.T_el = 6000.;
    liq_iron.Cv_el = 2.8;
    liq_iron.theta = 6000.;
    liq_iron.xi_0 = 15.786; 
    liq_iron.n = 1;
    
    return liq_iron;
  }
}
