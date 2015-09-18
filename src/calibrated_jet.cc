#include "calibrated_jet.h"

using namespace std;
using namespace fastjet;


MOD::CalibratedJet::CalibratedJet(double px, double py, double pz, double energy, string algorithm, double JEC, double area, double neutral_hadron_fraction, double neutral_em_fraction, int number_of_constituents, double charged_hadron_fraction, int charged_multiplicity, double charged_em_fraction) : _pseudojet(PseudoJet(px, py, pz, energy)), _algorithm(algorithm), _JEC(JEC), _area(area), _neutral_hadron_fraction(neutral_hadron_fraction), _neutral_em_fraction(neutral_em_fraction), _number_of_constituents(number_of_constituents), _charged_hadron_fraction(charged_hadron_fraction), _charged_multiplicity(charged_multiplicity), _charged_em_fraction(charged_em_fraction), _jet_quality(static_cast<JetQualityLevels_t>(-1)) {

}

MOD::CalibratedJet::CalibratedJet(PseudoJet pseudojet, string algorithm, double JEC, double area, double neutral_hadron_fraction, double neutral_em_fraction, int number_of_constituents, double charged_hadron_fraction, int charged_multiplicity, double charged_em_fraction) : _pseudojet(pseudojet), _algorithm(algorithm), _JEC(JEC), _area(area), _neutral_hadron_fraction(neutral_hadron_fraction), _neutral_em_fraction(neutral_em_fraction), _number_of_constituents(number_of_constituents), _charged_hadron_fraction(charged_hadron_fraction), _charged_multiplicity(charged_multiplicity), _charged_em_fraction(charged_em_fraction), _jet_quality(static_cast<JetQualityLevels_t>(-1)) {

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

   _jet_quality = static_cast<JetQualityLevels_t>(-1);
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
        << setw(16) << fixed << setprecision(8) << _pseudojet.px()
        << setw(16) << fixed << setprecision(8) << _pseudojet.py()
        << setw(16) << fixed << setprecision(8) << _pseudojet.pz()
        << setw(16) << fixed << setprecision(8) << _pseudojet.E()
        << setw(16) << fixed << setprecision(8) << _JEC 
        << setw(16) << fixed << setprecision(8) << _area
        << setw(16) << fixed << setprecision(8) << _neutral_hadron_fraction   
        << setw(16) << fixed << setprecision(8) << _neutral_em_fraction   
        << setw(16) << fixed << setprecision(8) << _number_of_constituents   
        << setw(16) << fixed << setprecision(8) << _charged_hadron_fraction   
        << setw(16) << fixed << setprecision(8) << _charged_multiplicity   
        << setw(16) << fixed << setprecision(8) << _charged_em_fraction       
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
   ss << "# AK5" << "              px              py              pz          energy             jec            area    neu_had_frac     neu_em_frac     no_of_const   chrg_had_frac     chrg_multip    chrg_em_frac" << endl;
   return ss.str();
}

string MOD::CalibratedJet::algorithm() const {
  return _algorithm;
}

PseudoJet MOD::CalibratedJet::uncorrected_pseudojet() const {
  return _pseudojet;
}

PseudoJet MOD::CalibratedJet::corrected_pseudojet() const {
  return _JEC * _pseudojet;
}

const int MOD::CalibratedJet::jet_quality() {

  // First, check if jet_quality has already been calculated or not. 

  if (_jet_quality != -1) {
    
    // This is the first time we're determining jet quality. Calculate it.
    
    double cut_offs [3] = { 0.99, 0.95, 0.90 };
    
    bool pass = false;

    for (unsigned i = 0; i < 3; i++) {
      
      pass = ( number_of_constituents() > 1 )     &&
             ( neutral_hadron_fraction() < cut_offs[i] ) && 
             ( neutral_em_fraction() < cut_offs[i] )     &&
             ( 
                ( abs(pseudojet().eta()) >= 2.4 ) || 
                ( charged_em_fraction() < 0.99 && charged_hadron_fraction() > 0.00 && charged_multiplicity() > 0) ); 
      
      // UNDETERMINED = -1, FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3

      if (pass) {
        _jet_quality =  static_cast<JetQualityLevels_t>(i + 1);
        return _jet_quality;
      }
    } 

    // If the code reached here, this means it didn't pass any quality cut. 
    _jet_quality = static_cast<JetQualityLevels_t>(0);
    
    return _jet_quality;
  }

  return _jet_quality;

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
