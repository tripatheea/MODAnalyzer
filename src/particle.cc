#include "../interface/particle.h"

using namespace std;
using namespace fastjet;


MODParticle::MODParticle(double px, double py, double pz, double energy, double mass, int pdgId, string trigger_type) : _pseudojet(PseudoJet(px, py, pz, energy)), _mass(mass), _pdgId(pdgId), _trigger_type(trigger_type) {
}

MODParticle::MODParticle(istringstream & input_stream) {

   string tag;
   double px, py, pz, energy, mass;
   int pdgId;

   input_stream >> tag >> px >> py >> pz >> energy >> mass >> pdgId;

   _trigger_type = tag;
   _pseudojet = PseudoJet(px, py, pz, energy);
   _mass = mass;
   _pdgId =pdgId;
}

MODParticle::MODParticle() {}

PseudoJet MODParticle::pseudojet() const {
   return _pseudojet;
}

int MODParticle::pdgId() const {
   return _pdgId;
}

double MODParticle::mass() const {
   return _mass;
}

string MODParticle::make_string() const {
   stringstream ss;
   ss << _trigger_type
        << setw(21) << setprecision(5) << _pseudojet.px()
        << setw(17) << setprecision(5) << _pseudojet.py()
        << setw(18) << setprecision(5) << _pseudojet.pz()
        << setw(18) << setprecision(5) << _pseudojet.E()
        << setw(19) << setprecision(5) << _mass
        << setw(18) << noshowpos << _pdgId
        << endl;

   return ss.str();
}

string MODParticle::make_header_string() const {
   stringstream ss;
   ss << "#" << _trigger_type << "               px               py               pz               energy               mass               pdgId" << endl;
   return ss.str();
}

ostream& operator<< (ostream& os, const MODParticle& particle) {
   os << particle.make_string();
   return os;
}
