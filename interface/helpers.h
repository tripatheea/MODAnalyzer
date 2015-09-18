#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>
#include <ostream>
#include <stdexcept>
#include <cmath>


#include "fastjet/ClusterSequence.hh"


namespace MOD {
  
  std::vector<MOD::PFCandidate> select_charged(std::vector<PFCandidate> pfcandidates); 
  std::vector<fastjet::PseudoJet> convert_to_pseudojets(std::vector<MOD::PFCandidate> pfcandidates);
  std::vector<fastjet::PseudoJet> convert_to_pseudojets(std::vector<MOD::CalibratedJet> jets);

  std::vector<fastjet::PseudoJet> filter_by_pT(std::vector<fastjet::PseudoJet>, double);


}