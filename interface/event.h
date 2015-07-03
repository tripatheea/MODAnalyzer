#include <iostream>
#include <vector>
#include <unordered_map>
#include <exception>
#include <memory>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <algorithm>


#include "trigger.h"
#include "pfcandidate.h"
#include "calibrated_jet.h"
#include "condition.h"

namespace MOD {

      class Event {

         public:
            Event(int, int, int, double);
            Event();

            int event_number() const;
            int run_number() const;

            int version() const;
            std::pair<std::string, std::string> data_type() const;

            double hardest_corrected_pt() const;
            double hardest_uncorrected_pt() const;

            const std::vector<PFCandidate> & particles() const;
            const std::vector<Trigger> & triggers() const;
            const std::vector<fastjet::PseudoJet> & pseudojets() const;

            const std::vector<CalibratedJet> corrected_calibrated_jets() const;          
            const std::vector<fastjet::PseudoJet> corrected_calibrated_pseudojets() const;

            const std::vector<CalibratedJet> & uncorrected_calibrated_jets() const;          
            const std::vector<fastjet::PseudoJet> & uncorrected_calibrated_pseudojets() const;

            std::string make_string() const;
            std::string assigned_trigger_name() const;

            const Trigger trigger_by_name(std::string name) const;    
            const Trigger assigned_trigger() const;

            void add_conditions(std::istringstream & input_stream); 
            void add_particle(std::istringstream & input_stream);
            void add_calibrated_jet(std::istringstream & input_stream);
            void add_trigger(std::istringstream & input_stream);
            
            void set_event_number(int event_number);
            void set_run_number(int run_number);
            void set_version(int version);
            void set_data_type(std::string a, std::string b);
            
            

            bool read_event(std::istream & data_stream);
            bool assigned_trigger_fired() const;

            int assigned_trigger_prescale() const;

            MOD::CalibratedJet hardest_corrected_jet();
            MOD::CalibratedJet hardest_uncorrected_jet();

            MOD::PFCandidate hardest_pfcandidate();

            friend std::ostream& operator<< (std::ostream&, const Event&);
            
         private:
            int _run_number, _event_number;


            

            double _hardest_corrected_pt;
            double _hardest_uncorrected_pt;

            
            int _version;
            std::pair<std::string, std::string> _data_type;

            std::string _assigned_trigger_name;
            
            Trigger _assigned_trigger;

            Condition _condition;

            std::vector<Trigger> _triggers;
            std::vector<PFCandidate> _particles;
            std::vector<fastjet::PseudoJet> _pseudojets;

            std::vector<CalibratedJet> _corrected_calibrated_jets;
            std::vector<fastjet::PseudoJet> _corrected_calibrated_pseudojets;

            std::vector<CalibratedJet> _uncorrected_calibrated_jets;
            std::vector<fastjet::PseudoJet> _uncorrected_calibrated_pseudojets;

            void set_assigned_trigger();
            void set_hardest_pt();
            void establish_properties();

      };

}