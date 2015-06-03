#include "../interface/fractional_jet_multiplicity.h"

using namespace std;
using namespace fastjet;


FractionalJetMultiplicity::FractionalJetMultiplicity(double cone_radius, double pt_cut) : _cone_radius(cone_radius), _pt_cut(pt_cut) {

}

const double FractionalJetMultiplicity::calculate_n_tilde(vector<PseudoJet> pseudojets) const {

	double N_tilde_current = 0.00;

	for(int i = 0; i < pseudojets.size(); i++) {
		double pt_i = pseudojets[i].pt();
		double pt_iR = 0.00;
		
		for(int j = 0; j < pseudojets.size(); j++) {
			double pt_j = pseudojets[j].pt();
			double squared_distance = pseudojets[i].squared_distance(pseudojets[j]);			// squared_distance instead of delta_R to speed things up.

			if (_cone_radius * _cone_radius > squared_distance)					// heavisideStep
				pt_iR += pt_j;
		}

		if (pt_iR > _pt_cut)	{							// heavisideStep
			N_tilde_current += pt_i / pt_iR;
		}
	}

	return N_tilde_current;
}