#ifndef EVENT_H
#define EVENT_H


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
#include <typeinfo>

#include "helpers.h"

#include "trigger.h"
#include "pfcandidate.h"
#include "calibrated_jet.h"
#include "condition.h"



namespace MOD {

      class Event {

         public:
            Event(int, int, int, double);
            Event();
            Event(int run_number, int Event_number, int version, std::pair<std::string, std::string> data_type, MOD::Condition condition, std::vector<MOD::Trigger> triggers, std::vector<MOD::PFCandidate> particles, std::vector<fastjet::PseudoJet> pseudojets, std::vector<MOD::CalibratedJet> CMS_jets, std::vector<fastjet::PseudoJet> CMS_pseudojets);

            int event_number() const;
            int run_number() const;

            int version() const;
            std::pair<std::string, std::string> data_type() const;

            const std::vector<PFCandidate> & particles() const;

            const std::vector<MOD::CalibratedJet> & CMS_jets() const;
            const std::vector<fastjet::PseudoJet> & CMS_pseudojets() const;
            const std::vector<fastjet::PseudoJet> & fastjet_pseudojets() const;

            const std::vector<Trigger> & triggers() const;

            std::string make_string() const;
            std::string assigned_trigger_name() const;

            const Trigger trigger_by_name(std::string name) const;    
            const Trigger assigned_trigger() const;

            void add_condition(std::istringstream & input_stream); 
            void add_particle(std::istringstream & input_stream);
            void add_CMS_jet(std::istringstream & input_stream);
            void add_trigger(std::istringstream & input_stream);
            
            void set_event_number(int event_number);
            void set_run_number(int run_number);
            void set_version(int version);
            void set_data_type(std::string a, std::string b);
            
            const bool trigger_exists(std::string trigger_name) const;
            bool read_event(std::istream & data_stream);
            bool assigned_trigger_fired() const;

            int assigned_trigger_prescale() const;

            // MOD::PFCandidate hardest_pfcandidate();


            fastjet::PseudoJet closest_fastjet_jet_to_trigger_jet();
            std::vector<fastjet::PseudoJet> closest_fastjet_jet_to_trigger_jet_constituents();

            bool trigger_jet_is_matched() const;

            const MOD::CalibratedJet trigger_jet() const;


            friend std::ostream& operator<< (std::ostream&, const Event&);
            
         private:
            int _run_number, _event_number;

            
            int _version;
            std::pair<std::string, std::string> _data_type;

            std::string _assigned_trigger_name;
            
            Trigger _assigned_trigger;

            Condition _condition;

            std::vector<MOD::Trigger> _triggers;
            std::vector<MOD::PFCandidate> _particles;
            std::vector<fastjet::PseudoJet> _pseudojets;

            std::vector<MOD::CalibratedJet> _CMS_jets;
            std::vector<fastjet::PseudoJet> _CMS_pseudojets;

            std::vector<fastjet::PseudoJet> _fastjet_pseudojets;

            MOD::CalibratedJet _trigger_jet;

            fastjet::PseudoJet _closest_fastjet_jet_to_trigger_jet;
            std::vector<fastjet::PseudoJet> _closest_fastjet_jet_to_trigger_jet_constituents;

            bool _trigger_jet_is_matched;

            void set_assigned_trigger();
            void set_hardest_pt();
            void establish_properties();

            void set_trigger_jet();
            void set_trigger_jet_is_matched();
            void set_closest_fastjet_jet_to_trigger_jet();

            std::vector<MOD::CalibratedJet> apply_jet_energy_corrections(std::vector<MOD::CalibratedJet> jets) const;
            std::vector<MOD::CalibratedJet> apply_eta_cut(std::vector<MOD::CalibratedJet> jets, double eta_cut) const;

      };
}


#endif /* EVENT_H */