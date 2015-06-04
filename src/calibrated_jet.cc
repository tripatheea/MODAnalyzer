#include "../interface/calibrated_jet.h"

using namespace std;
using namespace fastjet;


MOD::CalibratedJet::CalibratedJet(double px, double py, double pz, double energy, string algorithm) : _pseudojet(PseudoJet(px, py, pz, energy)), _algorithm(algorithm) {
}


MOD::CalibratedJet::CalibratedJet(istringstream & input_stream) {

   string tag;
   double px, py, pz, energy;

   input_stream >> tag >> px >> py >> pz >> energy;

   _algorithm = tag;
   _pseudojet = PseudoJet(px, py, pz, energy);
}

MOD::CalibratedJet::CalibratedJet() {}

PseudoJet MOD::CalibratedJet::pseudojet() const {
   return _pseudojet;
}

string MOD::CalibratedJet::make_string() const {
   stringstream ss;
   ss << _algorithm
        << setw(21) << setprecision(5) << _pseudojet.px()
        << setw(17) << setprecision(5) << _pseudojet.py()
        << setw(18) << setprecision(5) << _pseudojet.pz()
        << setw(18) << setprecision(5) << _pseudojet.E()
        << endl;

   return ss.str();
}

string MOD::CalibratedJet::make_header_string() const {
   stringstream ss;
   ss << "#" << _algorithm << "               px               py               pz               energy" << endl;
   return ss.str();
}

namespace MOD {
   ostream& operator<< (ostream& os, const CalibratedJet& jet) {
      os << jet.make_string();
      return os;
   }
}
