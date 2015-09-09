#include "calibrated_jet.h"

using namespace std;
using namespace fastjet;


MOD::CalibratedJet::CalibratedJet(double px, double py, double pz, double energy, string algorithm, double JEC, double area, double neutral_hadron_fraction, double neutral_em_fraction, int number_of_constituents, double charged_hadron_fraction, int charged_multiplicity, double charged_em_fraction) : _pseudojet(PseudoJet(px, py, pz, energy)), _algorithm(algorithm), _JEC(JEC), _area(area), _neutral_hadron_fraction(neutral_hadron_fraction), _neutral_em_fraction(neutral_em_fraction), _number_of_constituents(number_of_constituents), _charged_hadron_fraction(charged_hadron_fraction), _charged_multiplicity(charged_multiplicity), _charged_em_fraction(charged_em_fraction) {

}

MOD::CalibratedJet::CalibratedJet(PseudoJet pseudojet, string algorithm, double JEC, double area, double neutral_hadron_fraction, double neutral_em_fraction, int number_of_constituents, double charged_hadron_fraction, int charged_multiplicity, double charged_em_fraction) : _pseudojet(pseudojet), _algorithm(algorithm), _JEC(JEC), _area(area), _neutral_hadron_fraction(neutral_hadron_fraction), _neutral_em_fraction(neutral_em_fraction), _number_of_constituents(number_of_constituents), _charged_hadron_fraction(charged_hadron_fraction), _charged_multiplicity(charged_multiplicity), _charged_em_fraction(charged_em_fraction) {

}

MOD::CalibratedJet::CalibratedJet(istringstream & input_stream) {

   string tag;
   double px, py, pz, energy, JEC, area;

   double neutral_hadron_fraction, neutral_em_fraction, charged_hadron_fraction, charged_em_fraction;
   int number_of_constituents, charged_multiplicity;

   input_stream >> tag >> px >> py >> pz >> energy >> JEC >> area >> neutral_hadron_fraction >> neutral_em_fraction >> number_of_constituents >> charged_hadron_fraction >> charged_multiplicity >> charged_em_fraction;
   _pseudojet = PseudoJet(px, py, pz, energy);
   _algorithm = tag;
   _JEC = JEC;

   _area = area;

   _neutral_hadron_fraction = neutral_hadron_fraction;
   _neutral_em_fraction = neutral_em_fraction;
   _number_of_constituents = number_of_constituents;
   _charged_hadron_fraction = charged_hadron_fraction;
   _charged_multiplicity = charged_multiplicity;
   _charged_em_fraction = charged_em_fraction;

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
        << setw(20) << fixed << setprecision(8) << _area
        << setw(20) << fixed << setprecision(8) << _neutral_hadron_fraction   
        << setw(20) << fixed << setprecision(8) << _neutral_em_fraction   
        << setw(20) << fixed << setprecision(8) << _number_of_constituents   
        << setw(20) << fixed << setprecision(8) << _charged_hadron_fraction   
        << setw(20) << fixed << setprecision(8) << _charged_multiplicity   
        << setw(20) << fixed << setprecision(8) << _charged_em_fraction       
        << endl;

   return ss.str();
}

double MOD::CalibratedJet::JEC() const {
  return _JEC;
}

double MOD::CalibratedJet::area() const {
  return _area;
}


double MOD::CalibratedJet::neutral_hadron_fraction() const {
  return _neutral_hadron_fraction;
}

double MOD::CalibratedJet::neutral_em_fraction() const {
  return _neutral_em_fraction;
}

int MOD::CalibratedJet::number_of_constituents() const {
  return _number_of_constituents;
}

double MOD::CalibratedJet::charged_hadron_fraction() const {
  return _charged_hadron_fraction;
}

int MOD::CalibratedJet::charged_multiplicity() const {
  return _charged_multiplicity;
}

double MOD::CalibratedJet::charged_em_fraction() const {
  return _charged_em_fraction;
}

string MOD::CalibratedJet::make_header_string() const {
   stringstream ss;
   ss << "# AK5" << "                  px                  py                  pz              energy                 jec                area neutral_hadron_frac     neutral_em_frac  no_of_constituents charged_hadron_frac      charged_multip     charged_em_frac" << endl;
   return ss.str();
}

string MOD::CalibratedJet::algorithm() const {
  return _algorithm;
}

MOD::CalibratedJet MOD::CalibratedJet::corrected_jet() {
  PseudoJet new_pseudojet = _pseudojet * _JEC;

  MOD::CalibratedJet corrected_jet = MOD::CalibratedJet(new_pseudojet, _algorithm, 1.00, _area, _neutral_hadron_fraction, _neutral_em_fraction, _number_of_constituents, _charged_hadron_fraction, _charged_multiplicity, _charged_em_fraction);
  return corrected_jet;
}

bool MOD::CalibratedJet::jet_quality_cut(string level) {

   bool pass = false;
   double cut_off = 0.99;

   if (level == "tight") {
      cut_off = 0.90;
   }
   else if (level == "medium") {
      cut_off = 0.95;
   }
   else if (level == "loose") {
    cut_off = 0.99;
   }
   else {
    throw std::runtime_error("Invalid level for jet quality cut! Only 'loose', 'medium' and 'tight' are accepted.");
   }

   pass = ( number_of_constituents() > 1 )     &&
          ( neutral_hadron_fraction() < cut_off ) && 
          ( neutral_em_fraction() < cut_off )     &&
          ( 
            ( abs(pseudojet().eta()) >= 2.4 ) || 
               ( charged_em_fraction() < 0.99 && charged_hadron_fraction() > 0.00 && charged_multiplicity() > 0) ); 

   return pass;

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
