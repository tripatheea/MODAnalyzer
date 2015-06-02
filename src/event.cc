#include <iostream>
#include <unordered_map>
#include <exception>
#include <fstream>
#include <memory>
#include <string>
#include <iomanip> 

#include "fastjet/ClusterSequence.hh"
#include "trigger.cc"
#include "particle.cc"


using namespace std;
using namespace fastjet;

class MODEvent {

	public:
		MODEvent(int, int);
		MODEvent();

		int size();
		
		int event_number();
		int run_number();

		double calculate_N_tilde(double R, double pt_cut);	// R, pt_cut. R is the cone radius.

		vector<PseudoJet> jets(JetDefinition jet_def, double pt_cut);	// JetDefinition, pt_cut (Fastjet)
		vector<MODParticle> particles();

		void add_particle(double px, double py, double pz, double energy, double mass, int pdgId);
		void add_trigger(string name, int prescale_1, int prescale_2, bool fired);	
		void write_to_file(string filename);	// Will append if file already exists.

		string assigned_trigger_name();
		MODTrigger trigger_by_name(string name);

		double hardest_pt();

		void print_particles();

		vector<PseudoJet> particles_four_vectors();

		void assign_event_number(int MODEvent_number);
		void assign_run_number(int run_number);
		void assign_particles_type(string particles_type);

		bool is_valid();

		vector<MODTrigger> triggers();

		string type_of_particles();
		
	private:
		int _run_number, _event_number;
		vector<MODParticle> _particles;
		vector<MODTrigger> _triggers;
		string _type_of_particles;
};

MODEvent::MODEvent(int run_number, int MODEvent_number) : _run_number(run_number), _event_number(MODEvent_number) {}

MODEvent::MODEvent() {}

int MODEvent::event_number() {
	return _event_number;
}

int MODEvent::run_number() {
	return _run_number;
}


string MODEvent::type_of_particles() {
	return _type_of_particles;
}

void MODEvent::assign_run_number(int run_number) {
	_run_number = run_number;
}

void MODEvent::assign_particles_type(string particles_type) {
	_type_of_particles = particles_type;
}

void MODEvent::assign_event_number(int event_number) {
	_event_number = event_number;
}

void MODEvent::print_particles() {
	
}

int MODEvent::size() {
	return particles().size();
}

vector<PseudoJet> MODEvent::particles_four_vectors() {
	vector<PseudoJet> four_vectors;
	vector<MODParticle> all_particles = particles();

	for (unsigned int i = 0; i < all_particles.size(); i++) {
		MODParticle current_particle = all_particles[i];
		vector<double> current_particle_four_vector = current_particle.four_vector();
		four_vectors.push_back(PseudoJet(current_particle_four_vector[0], current_particle_four_vector[1], current_particle_four_vector[2], current_particle_four_vector[3]));
	}

	return four_vectors;
}

double MODEvent::calculate_N_tilde(double R, double pt_cut) {
	vector<PseudoJet> particles = particles_four_vectors();

	double N_tilde_current_MODEvent = 0.00;

	for(int i = 0; i < particles.size(); i++) {
		double pt_i = particles[i].pt();
		double pt_iR = 0.00;
		
		for(int j = 0; j < particles.size(); j++) {
			double pt_j = particles[j].pt();
			double squared_distance = particles[i].squared_distance(particles[j]);			// squared_distance instead of delta_R to speed things up.

			if (R*R > squared_distance)					// heavisideStep
				pt_iR += pt_j;
		}

		if (pt_iR > pt_cut)	{							// heavisideStep
			N_tilde_current_MODEvent += pt_i / pt_iR;
		}
	}

	return N_tilde_current_MODEvent;
}

vector<PseudoJet> MODEvent::jets(JetDefinition jet_def, double pt_cut) {
	vector<PseudoJet> particles = particles_four_vectors();

	// Run the clustering, extract the jets using fastjet.
	ClusterSequence cs(particles, jet_def);
	vector<PseudoJet> clustered_jets = cs.inclusive_jets(pt_cut);

	return clustered_jets;
}

vector<MODParticle> MODEvent::particles() {
	return _particles;
}

void MODEvent::add_particle(double px, double py, double pz, double energy, double mass, int pdgId) {
	MODParticle particle_to_add = MODParticle(px, py, pz, energy, mass, pdgId);
	_particles.push_back(particle_to_add);
}

void MODEvent::add_trigger(string name, int prescale_1, int prescale_2, bool fired) {
	pair <int, int> prescales;

	prescales = make_pair(prescale_1, prescale_2);

	MODTrigger trigger_to_add = MODTrigger(name, prescales, fired);
	_triggers.push_back(trigger_to_add);
}

MODTrigger MODEvent::trigger_by_name(string name) {
	vector<MODTrigger> triggers = this->triggers();

	for(int i = 0; i < triggers.size(); i++) {
		MODTrigger current_trigger = triggers[i];

		if (current_trigger.name() == name) {
			return current_trigger;
		}
	}

	MODTrigger * empty_trigger = new MODTrigger();
	return * empty_trigger;
}

vector<MODTrigger> MODEvent::triggers() {
	return _triggers;
}

void MODEvent::write_to_file(string filename) {
	ofstream file_to_write;

	file_to_write.open( filename, ios::out | ios::app ); 

	cout << "Writing things done!" << endl;

	int MODEvent_number = this->_event_number;
	int run_number = this->_run_number;

	vector<MODParticle> particles = this->particles();

	file_to_write << "BeginEvent Run " << this->_run_number << " Event " << this->_event_number << endl;
	
	// First, write out all particles.

	file_to_write << "#" << this->type_of_particles() << "               px               py               pz               energy               mass               pdgId" << endl;


	for (int i = 0; i < particles.size(); i++) {
		MODParticle current_particle = particles[i];

		vector<double> four_vector = current_particle.four_vector();

		file_to_write << type_of_particles() 
					  << setw(21) << setprecision(8) << four_vector[0] 
					  << setw(17) << setprecision(8) << four_vector[1] 
					  << setw(18) << setprecision(8) << four_vector[2] 
					  << setw(18) << setprecision(8) << four_vector[3] 
					  << setw(19) << setprecision(5) << current_particle.mass() 
					  << setw(18) << noshowpos << current_particle.pdgId() 
					  << endl;
	}

	// Next, write out all triggers.

	file_to_write << "#Trig          Name          Prescale_1          Prescale_2          Fired?" << endl;

	vector<MODTrigger> triggers = this->triggers();
	for(int i = 0; i < triggers.size(); i++) {
		MODTrigger current_trigger = triggers[i];

		pair<int, int> prescales = current_trigger.prescales();
		file_to_write << "trig" 
					  << setw(16) << current_trigger.name() 
					  << setw(15) << prescales.first 
					  << setw(20) << prescales.second 
					  << setw(17) << current_trigger.fired() 
					  << endl;
	}

	file_to_write << "EndEvent" << endl;
}

double MODEvent::hardest_pt() {
	vector<PseudoJet> particles = this->particles_four_vectors();

	// Run the clustering, extract the jets using fastjet.
	JetDefinition jet_def(antikt_algorithm, 0.5);
	ClusterSequence cs(particles, jet_def);
	vector<PseudoJet> clustered_jets = cs.inclusive_jets(0.0);
	
	double hardest_pt = 0.0;
	for (unsigned int i = 0; i < clustered_jets.size(); i++) {
		if (hardest_pt < clustered_jets[i].pt())
			hardest_pt = clustered_jets[i].pt();
	}

	return hardest_pt;
}

string MODEvent::assigned_trigger_name() {

	// Find the hardest jet first.
	



	JetDefinition jet_def(antikt_algorithm, 0.5);
	ClusterSequence cs(particles_four_vectors(), jet_def);
	vector<PseudoJet> clustered_jets = cs.inclusive_jets(0.0);;


	
	double hardest_pt = 0.0;
	for (unsigned int i = 0; i < clustered_jets.size(); i++) {
		if (hardest_pt < clustered_jets[i].pt())
			hardest_pt = clustered_jets[i].pt();
	}



	
	// Next, lookup which trigger to use based on the pt value of the hardest jet.
	
	/*
	37-56 GeV => 6U
	56-84 GeV => 15U
	84-114 GeV => 30U
	114-153 GeV => 50U
	>153 GeV => 70U
	*/

	string trigger_to_use;
	if (hardest_pt > 153) {
		trigger_to_use = "HLT_Jet70U";
	}
	else if (hardest_pt > 114) {
		trigger_to_use = "HLT_Jet50U";
	}
	else if (hardest_pt > 84) {
		trigger_to_use = "HLT_Jet30U";
	}
	else if (hardest_pt > 56) {
		trigger_to_use = "HLT_Jet15U";
	}
	else if (hardest_pt > 37) {
		trigger_to_use = "HLT_L1Jet6U";
	}
	else {
		trigger_to_use = "HLT_MinBiasPixel_SingleTrack";
	}
	

	// Here, we just return the trigger that was supposed to fire, not caring whether it actually did or not.
	// A check on whether it actually fired or not will be done in the N_tilde.cc file itself.

	// vector<string> triggersThatMatter {"HLT_L1Jet6U", "HLT_L1Jet10U", "HLT_Jet15U", "HLT_Jet30U", "HLT_Jet50U", "HLT_Jet70U", "HLT_Jet100U"};

	return trigger_to_use;


	

	// // Next, just see if the trigger_to_use fired or not.



	// if (trigger_to_use.length() != 0) {

	// 	Trigger selected_trigger = this->trigger_by_name(trigger_to_use);

	// 	if (selected_trigger.fired())
	// 		return selected_trigger;
	// }

	// // No trigger was fired for this MODEvent.
	// Trigger * empty_trigger = new Trigger();
	// return * empty_trigger;
}