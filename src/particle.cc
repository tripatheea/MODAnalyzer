#include "../interface/particle.h"

using namespace std;
using namespace fastjet;


MODParticle::MODParticle(double px, double py, double pz, double energy, double mass, int pdgId, string trigger_type) : _pseudojet(PseudoJet(px, py, pz, energy)), _mass(mass), _pdgId(pdgId), _trigger_type(trigger_type) {
}

MODParticle::MODParticle(string input_string) {
	vector<string> components = split(input_string);

	_trigger_type = components[0];
	_pseudojet = PseudoJet(stod(components[1]), stod(components[2]), stod(components[3]), stod(components[4]));
	_mass = stod(components[5]);
	_pdgId = stoi(components[6]);
}

MODParticle::MODParticle() {}

const PseudoJet MODParticle::pseudojet() const {
	return _pseudojet;
}

const int MODParticle::pdgId() const {
	return _pdgId;
}

const double MODParticle::mass() const {
	return _mass;
}

const string MODParticle::make_string() const {
	stringstream ss;

	ss << _trigger_type
		  << setw(21) << setprecision(8) << _pseudojet.px()
		  << setw(17) << setprecision(8) << _pseudojet.py()
		  << setw(18) << setprecision(8) << _pseudojet.pz()
		  << setw(18) << setprecision(8) << _pseudojet.E()
		  << setw(19) << setprecision(5) << _mass
		  << setw(18) << noshowpos << _pdgId
		  << endl;

	return ss.str();
}

const string MODParticle::header() const {
	stringstream ss;
	ss << "#" << _trigger_type << "               px               py               pz               energy               mass               pdgId" << endl;
	return ss.str();
}

vector<string> MODParticle::split(string const &input) { 
	istringstream buffer(input);
	vector<string> ret((istream_iterator<string>(buffer)), istream_iterator<string>());
	return ret;
}