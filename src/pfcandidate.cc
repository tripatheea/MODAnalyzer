#include "../interface/pfcandidate.h"

using namespace std;
using namespace fastjet;


MOD::PFCandidate::PFCandidate(double px, double py, double pz, double energy, double mass, int pdgId, string trigger_type) : _pdgId(pdgId), _trigger_type(trigger_type) {
   double recalc_energy = sqrt(px*px + py*py + pz*pz + mass*mass);

   if ( abs(recalc_energy - energy) > pow(10, -4)) {
      throw runtime_error("Recalculated energy (using 3-momentum nad mass) does not match give energy value.");
   }

   _pseudojet = PseudoJet(px, py, pz, recalc_energy);
}

MOD::PFCandidate::PFCandidate(istringstream & input_stream) {

   string tag;
   double px, py, pz, energy, mass;
   int pdgId;

   input_stream >> tag >> px >> py >> pz >> energy >> mass >> pdgId;
   
   double recalc_energy = sqrt(px*px + py*py + pz*pz + mass*mass);

   if ( abs(recalc_energy - energy) > pow(10, -4)) {
      throw runtime_error("Recalculated energy (using 3-momentum nad mass) does not match give energy value.");
   }

   _trigger_type = tag;
   _pseudojet = PseudoJet(px, py, pz, recalc_energy);
   _pdgId = pdgId;
}

MOD::PFCandidate::PFCandidate() {}

PseudoJet MOD::PFCandidate::pseudojet() const {
   return _pseudojet;
}

int MOD::PFCandidate::pdgId() const {
   return _pdgId;
}

double MOD::PFCandidate::mass() const {
   return _pseudojet.m();
}

string MOD::PFCandidate::make_string() const {
   stringstream ss;
   ss << _trigger_type
        << setw(21) << setprecision(5) << _pseudojet.px()
        << setw(17) << setprecision(5) << _pseudojet.py()
        << setw(18) << setprecision(5) << _pseudojet.pz()
        << setw(18) << setprecision(5) << _pseudojet.E()
        << setw(19) << setprecision(5) << mass()
        << setw(18) << noshowpos << _pdgId
        << endl;

   return ss.str();
}

string MOD::PFCandidate::make_header_string() const {
   stringstream ss;
   ss << "#" << _trigger_type << "               px               py               pz               energy               mass               pdgId" << endl;
   return ss.str();
}

namespace MOD {
   ostream& operator<< (ostream& os, const PFCandidate& particle) {
      os << particle.make_string();
      return os;
   }
}
