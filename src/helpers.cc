#include "helpers.h"

using namespace std;
using namespace fastjet;





std::vector<fastjet::PseudoJet> MOD::filter_charged(std::vector<fastjet::PseudoJet> particles) {
  vector<PFCandidate> filtered;

  for (unsigned i = 0; i < particles.size(); i++) {
    if ( (abs(particles[i].user_info<MOD::InfoPFC>().pdgId() == 211) || (abs(particles[i].user_info<MOD::InfoPFC>().pdgId() == 11) || (abs(particles[i].user_info<MOD::InfoPFC>().pdgId() == 13) ) {
      filtered.push_back(particles[i]);
    }
  }

  return filtered;
}

