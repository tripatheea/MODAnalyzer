#include "calibrated_jet.h"

using namespace std;
using namespace fastjet;


MOD::CalibratedJet::CalibratedJet(double px, double py, double pz, double energy, double mass, string algorithm, double JEC) : _pseudojet(PseudoJet(px, py, pz, energy)), _mass(mass), _algorithm(algorithm), _JEC(JEC) {
}


MOD::CalibratedJet::CalibratedJet(istringstream & input_stream) {

   string tag;
   double px, py, pz, energy, mass, JEC;

   input_stream >> tag >> px >> py >> pz >> energy >> mass >> JEC;

   _pseudojet = PseudoJet(px, py, pz, energy);
   _mass = mass;
   _algorithm = tag;
   _JEC = JEC;
}

MOD::CalibratedJet::CalibratedJet() {}

PseudoJet MOD::CalibratedJet::pseudojet() const {
   return _pseudojet;
}

string MOD::CalibratedJet::make_string() const {
   stringstream ss;
   ss << "  " << _algorithm
        << setw(14) << fixed << setprecision(8) << _pseudojet.px()
        << setw(14) << fixed << setprecision(8) << _pseudojet.py()
        << setw(14) << fixed << setprecision(8) << _pseudojet.pz()
        << setw(14) << fixed << setprecision(8) << _pseudojet.E()
        << setw(14) << fixed << setprecision(8) << _mass
        << setw(14) << fixed << setprecision(8) << _JEC
        << endl;

   return ss.str();
}

double MOD::CalibratedJet::JEC() const {
  return _JEC;
}

double MOD::CalibratedJet::mass() const {
  return _mass;
}

string MOD::CalibratedJet::make_header_string() const {
   stringstream ss;
   ss << "# " << _algorithm << "            px            py            pz        energy          mass           jec" << endl;
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
