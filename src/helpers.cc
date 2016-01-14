#include "helpers.h"


using namespace std;
using namespace fastjet;
// using namespace Pythia8;




std::vector<fastjet::PseudoJet> MOD::filter_charged(std::vector<fastjet::PseudoJet> particles) {
  vector<PseudoJet> filtered;


  Pythia8::Pythia pythia;

  for (unsigned i = 0; i < particles.size(); i++) {
    // if ( (abs(particles[i].user_info<MOD::InfoPFC>().pdgId()) == 211) || (abs(particles[i].user_info<MOD::InfoPFC>().pdgId()) == 11) || (abs(particles[i].user_info<MOD::InfoPFC>().pdgId()) == 13) ) {
    
    int pdgId = particles[i].user_info<MOD::InfoPFC>().pdgId();
    double charge = pythia.particleData.charge(pdgId);
    
    if ( charge != 0 ) {
      filtered.push_back(particles[i]);
    }


  }

  return filtered;
}

