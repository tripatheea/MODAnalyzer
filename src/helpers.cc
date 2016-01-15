#include "helpers.h"


using namespace std;
using namespace fastjet;




std::vector<fastjet::PseudoJet> MOD::filter_charged(std::vector<fastjet::PseudoJet> particles) {
  vector<PseudoJet> filtered;

  vector<int> pdg_ids {-11, -12, -13, -14, -16, -211, -2112, -2212, -321, 11, 12, 13, 130, 14, 16, 211, 2112, 22, 2212, 321, 32122};
  vector<int> charged_pdg_ids {11, 13, 15, 211, 321};

  for (unsigned i = 0; i < particles.size(); i++) {
    // if ( (abs(particles[i].user_info<MOD::InfoPFC>().pdgId()) == 211) || (abs(particles[i].user_info<MOD::InfoPFC>().pdgId()) == 11) || (abs(particles[i].user_info<MOD::InfoPFC>().pdgId()) == 13) ) {
    
    int pdgId = particles[i].user_info<MOD::InfoPFC>().pdgId();
    
    std::vector<int>::iterator it = find (charged_pdg_ids.begin(), charged_pdg_ids.end(), pdgId);

	if (it != charged_pdg_ids.end())
		filtered.push_back(particles[i]);

  }

  return filtered;
}

