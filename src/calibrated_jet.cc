#include "calibrated_jet.h"

using namespace std;
using namespace fastjet;


MOD::CalibratedJet::CalibratedJet(double px, double py, double pz, double energy, string algorithm, double JEC, double JEC_uncertainty, double area) : _pseudojet(PseudoJet(px, py, pz, energy)), _algorithm(algorithm), _JEC(JEC), _JEC_uncertainty(JEC_uncertainty), _area(area) {

}

MOD::CalibratedJet::CalibratedJet(PseudoJet pseudojet, string algorithm, double JEC, double JEC_uncertainty, double area) : _pseudojet(pseudojet), _algorithm(algorithm), _JEC(JEC), _JEC_uncertainty(JEC_uncertainty), _area(area) {

}

MOD::CalibratedJet::CalibratedJet(istringstream & input_stream) {

   string tag;
   double px, py, pz, energy, JEC, JEC_uncertainty, area;

   input_stream >> tag >> px >> py >> pz >> energy >> JEC >> JEC_uncertainty >> area;

   _pseudojet = PseudoJet(px, py, pz, energy);
   _algorithm = tag;
   _JEC = JEC;
   _JEC_uncertainty = JEC_uncertainty;
   _area = area;
}

MOD::CalibratedJet::CalibratedJet() {
}

bool MOD::CalibratedJet::is_valid() const {
   return ( ! _algorithm.empty());
}


PseudoJet MOD::CalibratedJet::pseudojet() const {
   return _pseudojet;
}

string MOD::CalibratedJet::make_string() const {
   stringstream ss;
   ss << "  " << _algorithm
        << setw(20) << fixed << setprecision(8) << _pseudojet.px()
        << setw(20) << fixed << setprecision(8) << _pseudojet.py()
        << setw(20) << fixed << setprecision(8) << _pseudojet.pz()
        << setw(20) << fixed << setprecision(8) << _pseudojet.E()
        << setw(20) << fixed << setprecision(8) << _JEC
        << setw(20) << fixed << setprecision(8) << _JEC_uncertainty        
        << setw(20) << fixed << setprecision(8) << _area
        << endl;

   return ss.str();
}

double MOD::CalibratedJet::JEC() const {
  return _JEC;
}

double MOD::CalibratedJet::JEC_uncertainty() const {
  return _JEC_uncertainty;
}

double MOD::CalibratedJet::area() const {
  return _area;
}

string MOD::CalibratedJet::make_header_string() const {
   stringstream ss;
   ss << "# AK5" << "                  px                  py                  pz              energy                 jec    jec_uncertainity                area" << endl;
   return ss.str();
}

string MOD::CalibratedJet::algorithm() const {
  return _algorithm;
}

MOD::CalibratedJet MOD::CalibratedJet::corrected_jet() {
  PseudoJet new_pseudojet = _pseudojet * _JEC;

  MOD::CalibratedJet corrected_jet = MOD::CalibratedJet(new_pseudojet, _algorithm, 1.00, _area, 0.0);
  return corrected_jet;
}


bool MOD::CalibratedJet::operator < (const MOD::CalibratedJet& j1) const {
  if (pseudojet().pt() > j1.pseudojet().pt())
    return true;
  return false;
}

namespace MOD {
   ostream& operator<< (ostream& os, const CalibratedJet& jet) {
      os << jet.make_string();
      return os;
   }

}
