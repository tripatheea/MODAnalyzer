#ifndef PFCANDIDATE_H
#define PFCANDIDATE_H

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
         PFCandidate(double px, double py, double pz, double energy, int pdgId);
         PFCandidate(std::istringstream & input_stream);
         PFCandidate(fastjet::PseudoJet particle);
         PFCandidate();

         fastjet::PseudoJet pseudojet() const;
         int pdgId() const;
         std::string make_string() const;
         std::string make_header_string() const;

         friend std::ostream& operator<< (std::ostream&, const PFCandidate&);

         bool operator<(const PFCandidate& rhs) const;

      private:
         fastjet::PseudoJet _pseudojet;
         int _pdgId; 
   };
}



#endif /* PFCANDIDATE_H */