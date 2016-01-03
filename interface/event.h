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
#include "condition.h"


#include "info_pfc.h"
#include "info_calibrated_jet.h"



#include "fastjet/Selector.hh"

namespace MOD {

      class Event {

         public:
            enum data_source_t { EXPERIMENT = 0, MC_TRUTH = 1, MC_RECO = 2, PRISTINE = 3 };


            Event(int, int, int, double);
            Event();

            const int event_number() const;
            const int run_number() const;

            const int version() const;
            const data_source_t data_source() const;
            const int weight() const;
            
            const std::pair<std::string, std::string> data_type() const;

            const std::vector<fastjet::PseudoJet> & particles() const;
            const std::vector<fastjet::PseudoJet> & cms_jets() const;
            const std::vector<fastjet::PseudoJet> & jets() const;

            const std::vector<Trigger> & triggers() const;

            std::string make_string() const;
            std::string assigned_trigger_name() const;

            const Trigger trigger_by_name(std::string name) const;    
            const Trigger trigger_by_short_name(std::string name) const;    
            const Trigger assigned_trigger() const;

            bool assigned_trigger_fired() const;
            int assigned_trigger_prescale() const;

            const fastjet::PseudoJet & closest_fastjet_jet_to_trigger_jet();
            const fastjet::PseudoJet & hardest_jet() const;
            const fastjet::PseudoJet & trigger_jet() const;

            bool trigger_jet_is_matched() const;
            
            void add_condition(std::istringstream & input_stream); 
            void add_particle(std::istringstream & input_stream);
            void add_cms_jet(std::istringstream & input_stream);
            void add_trigger(std::istringstream & input_stream);

            
            
            void convert_to_pristine();

            const bool trigger_exists(std::string trigger_name) const;  // You can give a trigger's full name or short name here.
            bool read_event(std::istream & data_stream);

            friend std::ostream& operator<< (std::ostream&, const Event&);
            
         private:
            int _run_number, _event_number, _version, _weight = 1;

            std::pair<std::string, std::string> _data_type;
            Condition _condition;

            std::vector<MOD::Trigger> _triggers;
            std::vector<fastjet::PseudoJet> _particles;
            std::vector<fastjet::PseudoJet> _cms_jets;
            std::vector<fastjet::PseudoJet> _jets;

            std::string _assigned_trigger_name;
            Trigger _assigned_trigger;
            
            fastjet::PseudoJet _trigger_jet;
            fastjet::PseudoJet _closest_fastjet_jet_to_trigger_jet;
            fastjet::PseudoJet _hardest_jet;

            bool _trigger_jet_is_matched;

            
            data_source_t _data_source;


            
            void set_assigned_trigger();
            void set_hardest_jet();
            void establish_properties();

            void set_trigger_jet();
            void set_trigger_jet_is_matched();
            void set_closest_fastjet_jet_to_trigger_jet();

            void set_event_number(int event_number);
            void set_run_number(int run_number);
            void set_version(int version);
            void set_data_type(std::string a, std::string b);
            void set_data_source(int data_source);

            std::vector<fastjet::PseudoJet> apply_jet_energy_corrections(std::vector<fastjet::PseudoJet> jets) const;

            const std::string stringify_jet(fastjet::PseudoJet jet) const;
            const std::string stringify_pfc(fastjet::PseudoJet particle) const;
      };
}


#endif /* EVENT_H */