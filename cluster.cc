#include "cluster.h"

using namespace std;
using namespace fastjet;


MODCluster::MODCluster(JetDefinition jet_definition, int pt_cut) : _jet_definition(jet_definition), _pt_cut(pt_cut) {

}

vector<PseudoJet> MODCluster::calculate_jets(MODEvent * event) {
	
	// Run the clustering, extract the jets using fastjet.

	ClusterSequence cs(event->particles_four_vectors(), _jet_definition);
	return cs.inclusive_jets(_pt_cut);
}