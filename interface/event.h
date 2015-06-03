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

		const int event_number() const;
		const int run_number() const;

		const double hardest_pt() const;

		const std::vector<MODParticle> & particles() const;
		const std::vector<MODTrigger> & triggers() const;
		const std::vector<fastjet::PseudoJet> & pseudojets() const;

		const std::string make_string();
		const std::string assigned_trigger_name() const;

		const MODTrigger & trigger_by_name(std::string name) const;		

		void add_particle(std::string input_string);
		void add_trigger(std::string input_string);	
		void set_event_number(int MODEvent_number);
		void set_run_number(int run_number);
		void set_particles_trigger_type(std::string trigger_type);
		
	private:
		int _run_number, _event_number;
				
		std::string _trigger_type;

		std::vector<MODParticle> _particles;
		std::vector<MODTrigger> _triggers;
		std::vector<PseudoJet> _pseudojets;

};