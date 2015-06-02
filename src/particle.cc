#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>

#include "fastjet/ClusterSequence.hh"

using namespace std;
using namespace fastjet;

class MODParticle {

	public:
		MODParticle(double px, double py, double pz, double energy, double mass, int pdgId, string trigger_type);
		MODParticle(string input_string);
		MODParticle();

		PseudoJet four_vector();
		int pdgId();
		double mass();
		string make_string();
		string header();

	private:
		string _trigger_type;
		PseudoJet _four_vector;
		double _mass;
		int _pdgId;	

		vector<string> split(string const &input);
};

MODParticle::MODParticle(double px, double py, double pz, double energy, double mass, int pdgId, string trigger_type) : _four_vector(PseudoJet(px, py, pz, energy)), _mass(mass), _pdgId(pdgId), _trigger_type(trigger_type) {
}

MODParticle::MODParticle(string input_string) {
	vector<string> components = this->split(input_string);

	_trigger_type = components[0];
	_four_vector = PseudoJet(stod(components[1]), stod(components[2]), stod(components[3]), stod(components[4]));
	_mass = stod(components[5]);
	_pdgId = stoi(components[6]);
}

MODParticle::MODParticle() {}

PseudoJet MODParticle::four_vector() {
	return _four_vector;
}

int MODParticle::pdgId() {
	return _pdgId;
}

double MODParticle::mass() {
	return _mass;
}

string MODParticle::make_string() {
	stringstream ss;

	ss << _trigger_type
		  << setw(21) << setprecision(8) << _four_vector.px()
		  << setw(17) << setprecision(8) << _four_vector.py()
		  << setw(18) << setprecision(8) << _four_vector.pz()
		  << setw(18) << setprecision(8) << _four_vector.E()
		  << setw(19) << setprecision(5) << _mass
		  << setw(18) << noshowpos << _pdgId
		  << endl;

	return ss.str();
}

string MODParticle::header() {
	stringstream ss;
	ss << "#" << _trigger_type << "               px               py               pz               energy               mass               pdgId" << endl;
	return ss.str();
}

vector<string> MODParticle::split(string const &input) { 
	istringstream buffer(input);
	vector<string> ret((istream_iterator<string>(buffer)), istream_iterator<string>());
	return ret;
}