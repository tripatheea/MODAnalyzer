#ifndef PDPFCANDIDATE_H
#define PDPFCANDIDATE_H

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
   class PDPFCandidate {

      public:
         PDPFCandidate(double px, double py, double pz, double energy, int pdgId);
         PDPFCandidate(std::istringstream & input_stream);
         PDPFCandidate();

         fastjet::PseudoJet pseudojet() const;
         int pdgId() const;
         std::string make_string() const;
         std::string make_header_string() const;

         friend std::ostream& operator<< (std::ostream&, const PDPFCandidate&);

         bool operator<(const PDPFCandidate& rhs) const;

      private:
         fastjet::PseudoJet _pseudojet;
         int _pdgId; 
   };
}



#endif /* PDPFCANDIDATE_H */