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

      int event_number() const;
      int run_number() const;

      double trigger_hardest_pt() const;

      const std::vector<MODParticle> & particles() const;
      const std::vector<MODTrigger> & triggers() const;
      const std::vector<fastjet::PseudoJet> & pseudojets() const;

      std::string make_string() const;
      std::string assigned_trigger_name() const;

      const MODTrigger & trigger_by_name(std::string name) const;    

      void add_particle(std::istringstream & input_stream);
      void add_trigger(std::istringstream & input_stream);  
      void set_event_number(int MODEvent_number);
      void set_run_number(int run_number);
      void set_particles_trigger_type(std::string trigger_type);

      bool read_event(ifstream & data_file);

      const MODTrigger & assigned_trigger() const;
      bool assigned_trigger_fired() const;
      int assigned_trigger_prescale() const;

      friend std::ostream& operator<< (std::ostream&, const MODEvent&);
      
   private:
      int _run_number, _event_number;

      double _trigger_hardest_pt;
      MODTrigger _assigned_trigger;
            
      std::string _trigger_type;

      std::vector<MODParticle> _particles;
      std::vector<MODTrigger> _triggers;
      std::vector<PseudoJet> _pseudojets;

      void set_assigned_trigger();
      void set_trigger_hardest_pt();
      void establish_properties();
};