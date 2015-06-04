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
   class PFCandidate {

      public:
         PFCandidate(double px, double py, double pz, double energy, double mass, int pdgId, std::string trigger_type);
         PFCandidate(std::istringstream & input_stream);
         PFCandidate();

         fastjet::PseudoJet pseudojet() const;
         int pdgId() const;
         double mass() const;
         std::string make_string() const;
         std::string make_header_string() const;

         friend std::ostream& operator<< (std::ostream&, const PFCandidate&);

      private:
         std::string _trigger_type;
         fastjet::PseudoJet _pseudojet;
         int _pdgId; 
   };
}