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






std::vector<fastjet::PseudoJet> MOD::convert_to_pseudojets(std::vector<MOD::CalibratedJet> jets) {
  vector<PseudoJet> pseudojets;

  for (unsigned i = 0; i < jets.size(); i++) {
    pseudojets.push_back(jets[i].uncorrected_pseudojet());
  }

  return pseudojets;
}


std::vector<fastjet::PseudoJet> MOD::convert_to_pseudojets(std::vector<MOD::MCCalibratedJet> mc_jets) {
  vector<PseudoJet> pseudojets;

  for (unsigned i = 0; i < mc_jets.size(); i++) {
    pseudojets.push_back(mc_jets[i].pseudojet());
  }

  return pseudojets;
}


std::vector<fastjet::PseudoJet> MOD::convert_to_pseudojets(std::vector<MOD::MCPFCandidate> mc_pfcandidates) {
  vector<PseudoJet> pseudojets;

  for (unsigned i = 0; i < mc_pfcandidates.size(); i++) {
    pseudojets.push_back(mc_pfcandidates[i].pseudojet());
  }

  return pseudojets;
}


std::vector<fastjet::PseudoJet> MOD::convert_to_pseudojets(std::vector<MOD::PDPFCandidate> pd_pfcandidates) {
  vector<PseudoJet> pseudojets;

  for (unsigned i = 0; i < pd_pfcandidates.size(); i++) {
    pseudojets.push_back(pd_pfcandidates[i].pseudojet());
  }

  return pseudojets;
}



std::vector<fastjet::PseudoJet> MOD::convert_to_pseudojets(std::vector<MOD::PDCalibratedJet> pd_jets) {
  vector<PseudoJet> pseudojets;

  for (unsigned i = 0; i < pd_jets.size(); i++) {
    pseudojets.push_back(pd_jets[i].pseudojet());
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

std::vector<MOD::PFCandidate> MOD::filter_by_pT(std::vector<MOD::PFCandidate> pfcandidates, double pT_lower_cut) {
  vector<PFCandidate> filtered;

  for (unsigned i = 0; i < pfcandidates.size(); i++) {
    if (pfcandidates[i].pseudojet().pt() > pT_lower_cut)
      filtered.push_back(pfcandidates[i]);
  }

  return filtered;
}



std::vector<fastjet::PseudoJet> MOD::filter_charged(std::vector<fastjet::PseudoJet> jet_constituents) {
  vector<PseudoJet> filtered_constituents;

  for (unsigned i = 0; i < jet_constituents.size(); i++) {
    if ( (abs(jet_constituents[i].user_index()) == 211) || (abs(jet_constituents[i].user_index()) == 11) || (abs(jet_constituents[i].user_index()) == 13) ) {
      filtered_constituents.push_back(jet_constituents[i]);
    }
  }

  return filtered_constituents;
}

std::vector<MOD::PFCandidate> MOD::filter_charged(std::vector<MOD::PFCandidate> pfcandidates) {
  vector<PFCandidate> filtered;

  for (unsigned i = 0; i < pfcandidates.size(); i++) {
    if ( (abs(pfcandidates[i].pseudojet().user_index()) == 211) || (abs(pfcandidates[i].pseudojet().user_index()) == 11) || (abs(pfcandidates[i].pseudojet().user_index()) == 13) ) {
      filtered.push_back(pfcandidates[i]);
    }
  }

  return filtered;
}

