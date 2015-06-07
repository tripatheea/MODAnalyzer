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
   class CalibratedJet {

      public:
         CalibratedJet(double px, double py, double pz, double energy, double mass, string algorithm);
         CalibratedJet(std::istringstream & input_stream);
         CalibratedJet();

         fastjet::PseudoJet pseudojet() const;
         
         std::string algorithm() const;
         std::string make_string() const;

         std::string make_header_string() const;

         friend std::ostream& operator<< (std::ostream&, const CalibratedJet&);

      private:
         std::string _algorithm;
         fastjet::PseudoJet _pseudojet;
         double _mass;
   };
}