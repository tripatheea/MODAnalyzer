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
#include "../src/pfcandidate.cc"
#include "../src/calibrated_jet.cc"

namespace MOD {

      class Event {

         public:
            Event(int, int);
            Event();

            int event_number() const;
            int run_number() const;

            double trigger_hardest_pt() const;

            const std::vector<PFCandidate> & particles() const;
            const std::vector<Trigger> & triggers() const;
            const std::vector<fastjet::PseudoJet> & pseudojets() const;
            
            const std::vector<fastjet::PseudoJet> & calibrated_jets_pseudojets() const;
            const std::vector<CalibratedJet> & calibrated_jets() const;

            std::string make_string() const;
            std::string assigned_trigger_name() const;

            const Trigger & trigger_by_name(std::string name) const;    

            void add_particle(std::istringstream & input_stream);
            void add_calibrated_jet(std::istringstream & input_stream);
            void add_trigger(std::istringstream & input_stream);

            void set_event_number(int event_number);
            void set_run_number(int run_number);
            void set_particles_trigger_type(std::string trigger_type);

            bool read_event(ifstream & data_file);

            const Trigger & assigned_trigger() const;
            bool assigned_trigger_fired() const;
            int assigned_trigger_prescale() const;

            friend std::ostream& operator<< (std::ostream&, const Event&);
            
         private:
            int _run_number, _event_number;

            double _trigger_hardest_pt;
            Trigger _assigned_trigger;
                  
            std::string _trigger_type;

            std::vector<Trigger> _triggers;
            std::vector<PFCandidate> _particles;
            std::vector<CalibratedJet> _calibrated_jets;
            std::vector<PseudoJet> _pseudojets;
            std::vector<PseudoJet> _calibrated_jets_pseudojets;

            void set_assigned_trigger();
            void set_trigger_hardest_pt();
            void establish_properties();
      };

}