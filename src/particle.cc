#include <iostream>
#include <vector>


using namespace std;

class MODParticle {

	public:
		MODParticle(double px, double py, double pz, double energy, double mass, int pdgId);
		MODParticle();

		vector<double> four_vector();
		int pdgId();
		double mass();

	private:
		vector<double> _four_vector;
		double _mass;
		int _pdgId;	
};

MODParticle::MODParticle(double px, double py, double pz, double energy, double mass, int pdgId) : _mass(mass), _pdgId(pdgId) {
	_four_vector.push_back(px);
	_four_vector.push_back(py);
	_four_vector.push_back(pz);
	_four_vector.push_back(energy);
}

MODParticle::MODParticle() {}

vector<double> MODParticle::four_vector() {
	return _four_vector;
}

int MODParticle::pdgId() {
	return _pdgId;
}

double MODParticle::mass() {
	return _mass;
}