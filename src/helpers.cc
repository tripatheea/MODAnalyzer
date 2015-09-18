#include "helpers.h"

using namespace std;
using namespace fastjet;




std::vector<MOD::PFCandidate> MOD::select_charged(std::vector<PFCandidate> pfcandidates) {

  vector<PFCandidate> charged_pfcandidates;

   for (unsigned i = 0; i < pfcandidates.size(); i++) {
      if ( (abs(pfcandidates[i].pdgId()) == 211) || (abs(pfcandidates[i].pdgId()) == 11) || (abs(pfcandidates[i].pdgId()) == 13) ) {
         charged_pfcandidates.push_back(pfcandidates[i]);
      }
   }

   return charged_pfcandidates;
}





std::vector<fastjet::PseudoJet> MOD::convert_to_pseudojets(std::vector<MOD::PFCandidate> pfcandidates) {
  vector<PseudoJet> pseudojets;

  for (unsigned i = 0; i < pfcandidates.size(); i++) {
    pseudojets.push_back(pfcandidates[i].pseudojet());
  }

  return pseudojets;
}




std::vector<fastjet::PseudoJet> MOD::convert_to_pseudojets(std::vector<MOD::CalibratedJet> jets) {
  vector<PseudoJet> pseudojets;

  for (unsigned i = 0; i < jets.size(); i++) {
    pseudojets.push_back(jets[i].uncorrected_pseudojet());
  }

  return pseudojets;
}


std::vector<fastjet::PseudoJet> MOD::filter_by_pT(std::vector<fastjet::PseudoJet> pseudojets, double pT_lower_cut) {
  vector<PseudoJet> filtered_pseudojets;
  
  for (unsigned i = 0; i < pseudojets.size(); i++) {
    if (pseudojets[i].pt() > pT_lower_cut)
      filtered_pseudojets.push_back(pseudojets[i]);
  }

  return filtered_pseudojets;
}


