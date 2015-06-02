#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 
#include <sstream>

using namespace std;
using namespace fastjet;


class MODNTilde {

	public:
		MODNTilde(double cone_radius, int pt_cut);
		double calculate_n_tilde(MODEvent * event);

	private:
		double _cone_radius;
		int _pt_cut;
};

MODNTilde::MODNTilde(double cone_radius, int pt_cut) : _cone_radius(cone_radius), _pt_cut(pt_cut) {

}

double MODNTilde::calculate_n_tilde(MODEvent * event) {

	vector<PseudoJet> particles = event->particles_four_vectors();

	double N_tilde_current_MODEvent = 0.00;

	for(int i = 0; i < particles.size(); i++) {
		double pt_i = particles[i].pt();
		double pt_iR = 0.00;
		
		for(int j = 0; j < particles.size(); j++) {
			double pt_j = particles[j].pt();
			double squared_distance = particles[i].squared_distance(particles[j]);			// squared_distance instead of delta_R to speed things up.

			if (_cone_radius * _cone_radius > squared_distance)					// heavisideStep
				pt_iR += pt_j;
		}

		if (pt_iR > _pt_cut)	{							// heavisideStep
			N_tilde_current_MODEvent += pt_i / pt_iR;
		}
	}

	return N_tilde_current_MODEvent;
}