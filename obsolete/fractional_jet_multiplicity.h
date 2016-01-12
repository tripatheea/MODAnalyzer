#ifndef FRACTIONAL_JET_MULTIPLICITY_H
#define FRACTIONAL_JET_MULTIPLICITY_H

#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>

#include "fastjet/ClusterSequence.hh"

namespace MOD {
   class FractionalJetMultiplicity {

      public:
         FractionalJetMultiplicity(double cone_radius, double pt_cut);
         double operator()(std::vector<fastjet::PseudoJet> pseudojets) const;

      private:
         const double _cone_radius;
         const double _pt_cut;

         double calculate_n_tilde(std::vector<fastjet::PseudoJet> pseudojets) const;
   };
}


#endif /* FRACTIONAL_JET_MULTIPLICITY_H */