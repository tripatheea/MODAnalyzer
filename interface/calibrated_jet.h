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
         CalibratedJet(double px, double py, double pz, double energy, std::string algorithm, double JEC, double JEC_uncertainty, double area);
         CalibratedJet(fastjet::PseudoJet pseudojet, std::string algorithm, double JEC, double JEC_uncertainty, double area);
         CalibratedJet(std::istringstream & input_stream);
         CalibratedJet();

         bool is_valid() const;

         fastjet::PseudoJet pseudojet() const;
         
         std::string algorithm() const;
         std::string make_string() const;

         std::string make_header_string() const;

         double JEC() const;
         double JEC_uncertainty() const;

         double area() const;

         MOD::CalibratedJet corrected_jet();

         friend std::ostream& operator<< (std::ostream&, const CalibratedJet&);

         bool operator<(const CalibratedJet& rhs) const;

      private:
         fastjet::PseudoJet _pseudojet;
         double _mass;
         std::string _algorithm;
         double _JEC;
         double _JEC_uncertainty;

         double _area;
   };
}