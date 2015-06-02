#include <iostream>
#include <vector>
#include <unordered_map>
#include <exception>
#include <fstream>
#include <memory>
#include <string>
#include <iomanip> 

#include "fastjet/ClusterSequence.hh"
#include "../src/trigger.cc"
#include "../src/particle.cc"


class MODEvent {

	public:
		MODEvent(int, int);
		MODEvent();

		int event_number();
		int run_number();

		std::vector<MODParticle> particles();
		std::vector<MODTrigger> triggers();
		std::vector<fastjet::PseudoJet> particles_four_vectors();

		void add_particle(std::string input_string);
		void add_trigger(std::string input_string);	
		void write_to_file(std::string filename);

		void set_event_number(int MODEvent_number);
		void set_run_number(int run_number);
		void set_particles_trigger_type(std::string trigger_type);

		double hardest_pt();

		std::string assigned_trigger_name();

		MODTrigger trigger_by_name(std::string name);		

	private:
		int _run_number, _event_number;
				
		std::string _trigger_type;

		std::vector<MODParticle> _particles;
		std::vector<MODTrigger> _triggers;

};