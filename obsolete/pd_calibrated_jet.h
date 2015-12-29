#ifndef PDCALIBRATED_JET_H
#define PDCALIBRATED_JET_H


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
   class PDCalibratedJet {

      public:
         PDCalibratedJet(double px, double py, double pz, double energy, std::string algorithm);
         PDCalibratedJet(fastjet::PseudoJet pseudojet, std::string algorithm);
         PDCalibratedJet(std::istringstream & input_stream);
         PDCalibratedJet();

         bool is_valid() const;
         
         std::string algorithm() const;
         std::string make_string() const;

         std::string make_header_string() const;

         fastjet::PseudoJet pseudojet() const;

         friend std::ostream& operator<< (std::ostream&, const PDCalibratedJet&);

         bool operator<(const PDCalibratedJet& rhs) const;

         bool operator==(const PDCalibratedJet& rhs) const;

      private:
         fastjet::PseudoJet _pseudojet;
         double _mass;
         std::string _algorithm;
   };
}



#endif /* PDCALIBRATED_JET_H */