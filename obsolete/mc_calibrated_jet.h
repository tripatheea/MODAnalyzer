#ifndef MCCALIBRATED_JET_H
#define MCCALIBRATED_JET_H


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
   class MCCalibratedJet {

      public:
         MCCalibratedJet(double px, double py, double pz, double energy, std::string algorithm);
         MCCalibratedJet(fastjet::PseudoJet pseudojet, std::string algorithm);
         MCCalibratedJet(std::istringstream & input_stream);
         MCCalibratedJet();

         bool is_valid() const;
         
         std::string algorithm() const;
         std::string make_string() const;

         std::string make_header_string() const;

         fastjet::PseudoJet pseudojet() const;

         friend std::ostream& operator<< (std::ostream&, const MCCalibratedJet&);

         bool operator<(const MCCalibratedJet& rhs) const;

         bool operator==(const MCCalibratedJet& rhs) const;

      private:
         fastjet::PseudoJet _pseudojet;
         double _mass;
         std::string _algorithm;
   };
}



#endif /* MCCALIBRATED_JET_H */