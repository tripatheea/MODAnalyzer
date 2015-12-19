#include "pd_calibrated_jet.h"

using namespace std;
using namespace fastjet;


MOD::PDCalibratedJet::PDCalibratedJet(double px, double py, double pz, double energy, string algorithm) {
  _pseudojet = PseudoJet(px, py, pz, energy);
  _algorithm = algorithm;
}

MOD::PDCalibratedJet::PDCalibratedJet(PseudoJet pseudojet, string algorithm) {
  _pseudojet = pseudojet;
  _algorithm = algorithm;
}

MOD::PDCalibratedJet::PDCalibratedJet(istringstream & input_stream) {

   string tag;
   double px, py, pz, energy;

   input_stream >> tag >> px >> py >> pz >> energy;
   _pseudojet = PseudoJet(px, py, pz, energy);
   _algorithm = tag;
}

MOD::PDCalibratedJet::PDCalibratedJet() {
}

bool MOD::PDCalibratedJet::is_valid() const {
   return ( ! _algorithm.empty());
}


string MOD::PDCalibratedJet::make_string() const {
   stringstream ss;
   ss << " PDAK5" << 
        << setw(16) << fixed << setprecision(8) << _pseudojet.px()
        << setw(16) << fixed << setprecision(8) << _pseudojet.py()
        << setw(16) << fixed << setprecision(8) << _pseudojet.pz()
        << setw(16) << fixed << setprecision(8) << _pseudojet.E()
        << endl;

   return ss.str();
}

string MOD::PDCalibratedJet::make_header_string() const {
   stringstream ss;
   ss << "# " << "PDAK5" << "" << "              px              py              pz          energy" << endl;
   return ss.str();
}

string MOD::PDCalibratedJet::algorithm() const {
  return _algorithm;
}


PseudoJet MOD::PDCalibratedJet::pseudojet() const {
  return _pseudojet;
}


bool MOD::PDCalibratedJet::operator < (const MOD::PDCalibratedJet& j1) const {
  if (pseudojet().pt() > j1.pseudojet().pt())
    return true;
  return false;
}

bool MOD::PDCalibratedJet::operator==(const PDCalibratedJet& rhs) const {
  return (_pseudojet == rhs._pseudojet) && (_algorithm == rhs._algorithm);
}

namespace MOD {
   ostream& operator<< (ostream& os, const PDCalibratedJet& jet) {
      os << jet.make_string();
      return os;
   }

}
