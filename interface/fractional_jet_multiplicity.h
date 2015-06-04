#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>

class FractionalJetMultiplicity {

   public:
      FractionalJetMultiplicity(double cone_radius, double pt_cut);
      double operator()(vector<PseudoJet> pseudojets) const;

   private:
      const double _cone_radius;
      const double _pt_cut;

      double calculate_n_tilde(vector<PseudoJet> pseudojets) const;
};