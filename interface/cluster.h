#include <vector>

class MODCluster {

	public:
		MODCluster(JetDefinition jet_definition, int pt_cut);
		vector<PseudoJet> calculate_jets(MODEvent * event);

	private:
		JetDefinition _jet_definition;
		int _pt_cut;
};
