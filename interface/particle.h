#include <vector>

class MODParticle {

	public:
		MODParticle(double px, double py, double pz, double energy, double mass, int pdgId, std::string trigger_type);
		MODParticle(std::string input_string);
		MODParticle();

		fastjet::PseudoJet four_vector();
		int pdgId();
		double mass();
		std::string make_string();
		std::string header();

	private:
		std::string _trigger_type;
		fastjet::PseudoJet _four_vector;
		double _mass;
		int _pdgId;	

		std::vector<std::string> split(std::string const &input);
};