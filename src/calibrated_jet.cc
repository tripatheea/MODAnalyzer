#include "../interface/calibrated_jet.h"

using namespace std;
using namespace fastjet;


MOD::CalibratedJet::CalibratedJet(double px, double py, double pz, double energy, double mass, string algorithm) : _pseudojet(PseudoJet(px, py, pz, energy)), _mass(mass), _algorithm(algorithm) {
}


MOD::CalibratedJet::CalibratedJet(istringstream & input_stream) {

   string tag;
   double px, py, pz, energy, mass;

   input_stream >> tag >> px >> py >> pz >> energy >> mass;

   _pseudojet = PseudoJet(px, py, pz, energy);
   _mass = mass;
   _algorithm = tag;
}

MOD::CalibratedJet::CalibratedJet() {}

PseudoJet MOD::CalibratedJet::pseudojet() const {
   return _pseudojet;
}

string MOD::CalibratedJet::make_string() const {
   stringstream ss;
   ss << "  " << _algorithm
        << setw(12) << fixed << setprecision(5) << _pseudojet.px()
        << setw(12) << fixed << setprecision(5) << _pseudojet.py()
        << setw(12) << fixed << setprecision(5) << _pseudojet.pz()
        << setw(12) << fixed << setprecision(5) << _pseudojet.E()
        << setw(12) << fixed << setprecision(5) << _mass
        << endl;

   return ss.str();
}

string MOD::CalibratedJet::make_header_string() const {
   stringstream ss;
   ss << "# " << _algorithm << "          px          py          pz      energy        mass" << endl;
   return ss.str();
}

string MOD::CalibratedJet::algorithm() const {
  return _algorithm;
}

namespace MOD {
   ostream& operator<< (ostream& os, const CalibratedJet& jet) {
      os << jet.make_string();
      return os;
   }
}
