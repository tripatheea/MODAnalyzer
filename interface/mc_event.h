#ifndef MCEVENT_H
#define MCEVENT_H


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

      class MCEvent {

         public:
            MCEvent(int, int, int, double);
            MCEvent();
            MCEvent(int run_number, int Event_number, int version, std::pair<std::string, std::string> data_type, MOD::Condition condition, std::vector<MOD::Trigger> triggers, std::vector<MOD::PFCandidate> particles, std::vector<fastjet::PseudoJet> pseudojets, std::vector<MOD::CalibratedJet> CMS_jets, std::vector<fastjet::PseudoJet> CMS_pseudojets);

            int event_number() const;
            int run_number() const;

            int version() const;
            std::pair<std::string, std::string> data_type() const;

            const std::vector<PFCandidate> & particles() const;

            const std::vector<fastjet::PseudoJet> & pseudojets() const;


            std::string make_string() const;

            

            void add_condition(std::istringstream & input_stream); 
            void add_particle(std::istringstream & input_stream);
            
            
            void set_event_number(int event_number);
            void set_run_number(int run_number);
            void set_version(int version);
            void set_data_type(std::string a, std::string b);
            
            
            bool read_event(std::istream & data_stream);
            
                        
            const fastjet::PseudoJet hardest_pseudojet() const;


            friend std::ostream& operator<< (std::ostream&, const MCEvent&);
            
         private:
            int _run_number, _event_number;

            
            int _version;
            std::pair<std::string, std::string> _data_type;

            

            Condition _condition;

            std::vector<MOD::PFCandidate> _particles;
            std::vector<fastjet::PseudoJet> _pseudojets;

            
            std::vector<fastjet::PseudoJet> _fastjet_pseudojets;



            void set_assigned_trigger();
            void set_hardest_pt();
            void establish_properties();

            std::vector<MOD::CalibratedJet> apply_eta_cut(std::vector<MOD::CalibratedJet> jets, double eta_cut) const;

      };
}


#endif /* MCEVENT_H */