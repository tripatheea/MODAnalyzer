#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>

class FractionalJetMultiplicity {

	public:
		FractionalJetMultiplicity(double cone_radius, double pt_cut);
		double calculate_n_tilde(vector<PseudoJet> pseudojets);

	private:
		double _cone_radius;
		double _pt_cut;
};