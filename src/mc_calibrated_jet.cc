#include "mc_calibrated_jet.h"

using namespace std;
using namespace fastjet;


MOD::MCCalibratedJet::MCCalibratedJet(double px, double py, double pz, double energy, string algorithm) {

}

MOD::MCCalibratedJet::MCCalibratedJet(PseudoJet pseudojet, string algorithm) {

}

MOD::MCCalibratedJet::MCCalibratedJet(istringstream & input_stream) {

   string tag;
   double px, py, pz, energy;

   input_stream >> tag >> px >> py >> pz >> energy;
   _pseudojet = PseudoJet(px, py, pz, energy);
   _algorithm = tag;
}

MOD::MCCalibratedJet::MCCalibratedJet() {
}

bool MOD::MCCalibratedJet::is_valid() const {
   return ( ! _algorithm.empty());
}


string MOD::MCCalibratedJet::make_string() const {
   stringstream ss;
   ss << "  " << _algorithm
        << setw(16) << fixed << setprecision(8) << _pseudojet.px()
        << setw(16) << fixed << setprecision(8) << _pseudojet.py()
        << setw(16) << fixed << setprecision(8) << _pseudojet.pz()
        << setw(16) << fixed << setprecision(8) << _pseudojet.E()
        << endl;

   return ss.str();
}

string MOD::MCCalibratedJet::make_header_string() const {
   stringstream ss;
   ss << "# " << _algorithm << "" << "              px              py              pz          energy" << endl;
   return ss.str();
}

string MOD::MCCalibratedJet::algorithm() const {
  return _algorithm;
}


PseudoJet MOD::MCCalibratedJet::pseudojet() const {
  return _pseudojet;
}


bool MOD::MCCalibratedJet::operator < (const MOD::MCCalibratedJet& j1) const {
  if (pseudojet().pt() > j1.pseudojet().pt())
    return true;
  return false;
}

bool MOD::MCCalibratedJet::operator==(const MCCalibratedJet& rhs) const {
  return (_pseudojet == rhs._pseudojet) && (_algorithm == rhs._algorithm);
}

namespace MOD {
   ostream& operator<< (ostream& os, const MCCalibratedJet& jet) {
      os << jet.make_string();
      return os;
   }

}
