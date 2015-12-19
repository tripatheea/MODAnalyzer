#ifndef HELPERS_H
#define HELPERS_H


#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>
#include <ostream>
#include <stdexcept>
#include <cmath>


#include "fastjet/ClusterSequence.hh"
#include "pfcandidate.h"
#include "calibrated_jet.h"

#include "mc_pfcandidate.h"
#include "mc_calibrated_jet.h"




namespace MOD {
  
  std::vector<MOD::PFCandidate> select_charged(std::vector<MOD::PFCandidate> pfcandidates); 
  
  std::vector<fastjet::PseudoJet> convert_to_pseudojets(std::vector<MOD::PFCandidate> pfcandidates);
  std::vector<fastjet::PseudoJet> convert_to_pseudojets(std::vector<MOD::CalibratedJet> jets);

  std::vector<fastjet::PseudoJet> convert_to_pseudojets(std::vector<MOD::MCPFCandidate> mc_pfcandidates);
  std::vector<fastjet::PseudoJet> convert_to_pseudojets(std::vector<MOD::MCCalibratedJet> mc_jets);

  std::vector<fastjet::PseudoJet> filter_by_pT(std::vector<fastjet::PseudoJet>, double pT_cut);
  std::vector<MOD::PFCandidate> filter_by_pT(std::vector<MOD::PFCandidate>, double pT_cut);

  std::vector<fastjet::PseudoJet> filter_charged(std::vector<fastjet::PseudoJet>);
  std::vector<MOD::PFCandidate> filter_charged(std::vector<MOD::PFCandidate>);

}


#endif /* HELPERS_H */