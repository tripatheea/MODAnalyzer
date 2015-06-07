#include <iostream>
#include <vector>
#include <unordered_map>
#include <exception>
#include <memory>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <limits>



#include "../interface/trigger.h"
#include "../interface/pfcandidate.h"
#include "../interface/calibrated_jet.h"

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
            const std::vector<CalibratedJet> & calibrated_jets_ak5() const;          
            const std::vector<fastjet::PseudoJet> & calibrated_pseudojets_ak5() const;
            const std::vector<CalibratedJet> & calibrated_jets_ak7() const;
            const std::vector<fastjet::PseudoJet> & calibrated_pseudojets_ak7() const;

            std::string make_string() const;
            std::string assigned_trigger_name() const;

            const Trigger & trigger_by_name(std::string name) const;    
            const Trigger & assigned_trigger() const;

            void add_particle(std::istringstream & input_stream);
            void add_calibrated_jet(std::istringstream & input_stream);
            void add_trigger(std::istringstream & input_stream);
            void set_event_number(int event_number);
            void set_run_number(int run_number);

            bool read_event(std::istream & data_stream);
            bool assigned_trigger_fired() const;

            int assigned_trigger_prescale() const;

            friend std::ostream& operator<< (std::ostream&, const Event&);
            
         private:
            int _run_number, _event_number;

            std::string _assigned_trigger_name;

            double _trigger_hardest_pt;
            Trigger _assigned_trigger;

            std::vector<Trigger> _triggers;
            std::vector<PFCandidate> _particles;
            std::vector<fastjet::PseudoJet> _pseudojets;

            std::vector<CalibratedJet> _calibrated_jets_ak5;
            std::vector<fastjet::PseudoJet> _calibrated_pseudojets_ak5;

            std::vector<CalibratedJet> _calibrated_jets_ak7;
            std::vector<fastjet::PseudoJet> _calibrated_pseudojets_ak7;

            void set_assigned_trigger();
            void set_trigger_hardest_pt();
            void establish_properties();
      };

}