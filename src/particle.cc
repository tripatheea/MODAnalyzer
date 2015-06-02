#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>



using namespace std;

class MODParticle {

	public:
		MODParticle(double px, double py, double pz, double energy, double mass, int pdgId, string trigger_type);
		MODParticle();

		vector<double> four_vector();
		int pdgId();
		double mass();
		string make_string();
		string header();

	private:
		string _trigger_type;
		double _px; double _py; double _pz; double _energy;
		double _mass;
		int _pdgId;	
};

MODParticle::MODParticle(double px, double py, double pz, double energy, double mass, int pdgId, string trigger_type) : _px(px), _py(py), _pz(pz), _energy(energy), _mass(mass), _pdgId(pdgId), _trigger_type(trigger_type) {
}

MODParticle::MODParticle() {}

vector<double> MODParticle::four_vector() {
	vector<double> four_vector = {_px, _py, _pz, _energy};
	return four_vector;
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
		  << setw(21) << setprecision(8) << _px
		  << setw(17) << setprecision(8) << _py
		  << setw(18) << setprecision(8) << _pz
		  << setw(18) << setprecision(8) << _energy
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