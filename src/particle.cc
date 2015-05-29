#include <iostream>
#include <vector>


using namespace std;

class Particle {

	public:
		Particle(double px, double py, double pz, double energy, double mass, int pdgId);
		Particle();

		vector<double> four_vector();
		int pdgId();
		double mass();

	private:
		vector<double> four_vector_;
		double mass_;
		int pdgId_;	
};

Particle::Particle(double px, double py, double pz, double energy, double mass, int pdgId) : mass_(mass), pdgId_(pdgId) {
	four_vector_.push_back(px);
	four_vector_.push_back(py);
	four_vector_.push_back(pz);
	four_vector_.push_back(energy);
}

Particle::Particle() {}

vector<double> Particle::four_vector() {
	return four_vector_;
}

int Particle::pdgId() {
	return Particle::pdgId_;
}

double Particle::mass() {
	return Particle::mass_;
}