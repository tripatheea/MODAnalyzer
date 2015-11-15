#ifndef MCMCPFCANDIDATE_H
#define MCPFCANDIDATE_H

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
   class MCPFCandidate {

      public:
         MCPFCandidate(double px, double py, double pz, double energy, int pdgId);
         MCPFCandidate(std::istringstream & input_stream);
         MCPFCandidate();

         fastjet::PseudoJet pseudojet() const;
         int pdgId() const;
         std::string make_string() const;
         std::string make_header_string() const;

         friend std::ostream& operator<< (std::ostream&, const MCPFCandidate&);

         bool operator<(const MCPFCandidate& rhs) const;

      private:
         fastjet::PseudoJet _pseudojet;
         int _pdgId; 
   };
}



#endif /* MCMCPFCANDIDATE_H */