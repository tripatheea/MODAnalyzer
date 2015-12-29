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
            Event(int, int, int, double);
            Event();

            int event_number() const;
            int run_number() const;

            int version() const;
            std::pair<std::string, std::string> data_type() const;

            const std::vector<fastjet::PseudoJet> & particles() const;
            const std::vector<fastjet::PseudoJet> & cms_jets() const;
            const std::vector<fastjet::PseudoJet> & jets() const;

            const std::vector<Trigger> & triggers() const;

            std::string make_string() const;
            std::string assigned_trigger_name() const;

            const Trigger trigger_by_name(std::string name) const;    
            const Trigger trigger_by_short_name(std::string name) const;    
            const Trigger assigned_trigger() const;

            void add_condition(std::istringstream & input_stream); 
            
            void add_particle(std::istringstream & input_stream);
            void add_cms_jet(std::istringstream & input_stream);
            
            void add_trigger(std::istringstream & input_stream);

            
            
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
            const fastjet::PseudoJet hardest_jet() const;


          
            bool trigger_jet_is_matched() const;

            const fastjet::PseudoJet trigger_jet() const;

            const int data_source() const;
            void set_data_source(int data_source);

      
            const int weight() const;
            
            void convert_to_pristine();

            friend std::ostream& operator<< (std::ostream&, const Event&);
            
         private:
            int _run_number, _event_number;

            
            int _version;
            std::pair<std::string, std::string> _data_type;

            std::string _assigned_trigger_name;
            
            Trigger _assigned_trigger;

            Condition _condition;


            
            std::vector<MOD::Trigger> _triggers;

            std::vector<fastjet::PseudoJet> _particles;

            std::vector<fastjet::PseudoJet> _cms_jets;
            std::vector<fastjet::PseudoJet> _jets;

            fastjet::PseudoJet _trigger_jet;
            fastjet::PseudoJet _closest_fastjet_jet_to_trigger_jet;
            fastjet::PseudoJet _hardest_jet;

            bool _trigger_jet_is_matched;

            
            
            
            void set_assigned_trigger();
            void set_hardest_jet();
            void establish_properties();

            void set_trigger_jet();
            void set_trigger_jet_is_matched();
            void set_closest_fastjet_jet_to_trigger_jet();


            void set_weight(int weight);

            std::vector<fastjet::PseudoJet> apply_jet_energy_corrections(std::vector<fastjet::PseudoJet> jets) const;
            
            enum data_source_t { EXPERIMENT = 0, MC_TRUTH = 1, MC_RECO = 2, PRISTINE = 3 };
            data_source_t _data_source;

            int _weight = 1;

            
      };
}


#endif /* EVENT_H */