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
// #include "mc_pfcandidate.h"
#include "calibrated_jet.h"
#include "mc_calibrated_jet.h"
#include "condition.h"



namespace MOD {

      class Event {

         public:
            Event(int, int, int, double);
            Event();
            // Event(int run_number, int Event_number, int version, std::pair<std::string, std::string> data_type, MOD::Condition condition, std::vector<MOD::Trigger> triggers, std::vector<MOD::PFCandidate> particles, std::vector<MOD::CalibratedJet> CMS_jets);

            int event_number() const;
            int run_number() const;

            int version() const;
            std::pair<std::string, std::string> data_type() const;

            const std::vector<PFCandidate> & particles() const;

            const std::vector<MOD::CalibratedJet> & CMS_jets() const;
            const std::vector<fastjet::PseudoJet> & fastjet_clustered_pseudojets() const;

            const std::vector<Trigger> & triggers() const;

            std::string make_string() const;
            std::string assigned_trigger_name() const;

            const Trigger trigger_by_name(std::string name) const;    
            const Trigger trigger_by_short_name(std::string name) const;    
            const Trigger assigned_trigger() const;

            void add_condition(std::istringstream & input_stream); 
            void add_particle(std::istringstream & input_stream);
            void add_CMS_jet(std::istringstream & input_stream);
            void add_trigger(std::istringstream & input_stream);

            void add_mc_truth_particle(std::istringstream & input_stream);
            void add_mc_reco_particle(std::istringstream & input_stream);
            
            void set_event_number(int event_number);
            void set_run_number(int run_number);
            void set_version(int version);
            void set_data_type(std::string a, std::string b);
            
            // You can give a trigger's full name or short name here.
            const bool trigger_exists(std::string trigger_name) const;
            bool read_event(std::istream & data_stream);
            bool assigned_trigger_fired() const;

            int assigned_trigger_prescale() const;



            fastjet::PseudoJet closest_fastjet_jet_to_trigger_jet();

            const MOD::MCCalibratedJet hardest_mc_truth_jet() const;
            const MOD::MCCalibratedJet hardest_mc_reco_jet() const;




            void set_hardest_truth_jet();
            void set_hardest_reco_jet();


            bool trigger_jet_is_matched() const;

            const MOD::CalibratedJet trigger_jet() const;

            const int data_source() const;
            void set_data_source(int data_source);

            const bool pristine_form() const;
            void set_pristine_form(bool pristine);

            const int prescale() const;
            void set_prescale(int prescale);


            friend std::ostream& operator<< (std::ostream&, const Event&);
            
         private:
            int _run_number, _event_number;

            
            int _version;
            std::pair<std::string, std::string> _data_type;

            std::string _assigned_trigger_name;
            
            Trigger _assigned_trigger;

            Condition _condition;


            // Things related to CMS Data.
            std::vector<MOD::Trigger> _triggers;
            std::vector<MOD::PFCandidate> _particles;

            std::vector<MOD::CalibratedJet> _CMS_jets;
            std::vector<fastjet::PseudoJet> _fastjet_clustered_pseudojets;

            MOD::CalibratedJet _trigger_jet;

            fastjet::PseudoJet _closest_fastjet_jet_to_trigger_jet;
            

            bool _trigger_jet_is_matched;


            // Things related to Monte Carlo.
            // Note that we don't directly read AK5 jets because in the Monte Carlo generated MOD files, we don't write out AK5 jets. 
            std::vector<MOD::MCPFCandidate> _mc_truth_particles;
            std::vector<MOD::MCPFCandidate> _mc_reco_particles;

            std::vector<MOD::MCCalibratedJet> _mc_truth_jets;
            std::vector<MOD::MCCalibratedJet> _mc_reco_jets;

            MOD::MCCalibratedJet _hardest_mc_truth_jet;
            MOD::MCCalibratedJet _hardest_mc_reco_jet;

            

            
            void set_assigned_trigger();
            void set_hardest_pt();
            void establish_properties();

            void set_trigger_jet();
            void set_trigger_jet_is_matched();
            void set_closest_fastjet_jet_to_trigger_jet();

            std::vector<MOD::CalibratedJet> apply_jet_energy_corrections(std::vector<MOD::CalibratedJet> jets) const;
            std::vector<MOD::CalibratedJet> apply_eta_cut(std::vector<MOD::CalibratedJet> jets, double eta_cut) const;

            enum data_source_t { EXPERIMENT = 0, MC_TRUTH = 1, MC_RECO = 2 };
            data_source_t _data_source;

            bool _pristine_form;

            int _prescale = 1;

            
      };
}


#endif /* EVENT_H */