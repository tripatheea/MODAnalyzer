#include "../interface/event.h"

using namespace std;
using namespace fastjet;

MODEvent::MODEvent(int run_number, int MODEvent_number) : _run_number(run_number), _event_number(MODEvent_number) {}

MODEvent::MODEvent() {}

int MODEvent::event_number() const {
   return _event_number;
}

int MODEvent::run_number() const {
   return _run_number;
}

void MODEvent::set_run_number(int run_number) {
   _run_number = run_number;
}

void MODEvent::set_particles_trigger_type(string trigger_type) {
   _trigger_type = trigger_type;
}

void MODEvent::set_event_number(int event_number) {
   _event_number = event_number;
}

const vector<PseudoJet> & MODEvent::pseudojets() const {
   return _pseudojets;
}

const vector<MODParticle> & MODEvent::particles() const {
   return _particles;
}

void MODEvent::add_particle(istringstream & input_stream) {
   MODParticle new_particle = MODParticle(input_stream);
   _particles.push_back(new_particle);
   _pseudojets.push_back(PseudoJet(new_particle.pseudojet()));
}

void MODEvent::add_trigger(istringstream & input_stream) {
   _triggers.push_back(MODTrigger(input_stream));
}

const MODTrigger & MODEvent::trigger_by_name(string name) const {
   for(int i = 0; i < triggers().size(); i++) {
      const MODTrigger& current_trigger = triggers().at(i);

      if (current_trigger.name() == name) {
         return current_trigger;
      }
   }

   MODTrigger * empty_trigger = new MODTrigger();
   return *empty_trigger;
}

const vector<MODTrigger> & MODEvent::triggers() const {
   return _triggers;
}

string MODEvent::make_string() const {
   stringstream file_to_write;
   
   file_to_write << "BeginEvent Run " << _run_number << " Event " << _event_number << endl;
   
   // First, write out all particles.

   file_to_write << _particles[0].make_header_string();
   for (int i = 0; i < _particles.size(); i++) {
      file_to_write << _particles[i];
   }

   // Next, write out all triggers.

   file_to_write << _triggers[0].make_header_string();
   for(int i = 0; i < _triggers.size(); i++) {
      file_to_write << _triggers[i];
   }

   file_to_write << "EndEvent" << endl;

   return file_to_write.str();
}

double MODEvent::trigger_hardest_pt() const {

   // Run the clustering, extract the jets using fastjet.
   JetDefinition jet_def(antikt_algorithm, 0.5);
   ClusterSequence cs(pseudojets(), jet_def);
   vector<PseudoJet> clustered_jets = cs.inclusive_jets(0.0);

   double hardest_pt = 0.0;
   for (unsigned int i = 0; i < clustered_jets.size(); i++) {
      if (hardest_pt < clustered_jets[i].pt()) {
         hardest_pt = clustered_jets[i].pt();
      }
   }

   return hardest_pt;
}

string MODEvent::assigned_trigger_name() const {

   double hardest_pt_value = trigger_hardest_pt();

   // Next, lookup which trigger to use based on the pt value of the hardest jet.

   string trigger_to_use;
   if (hardest_pt_value > 153) {
      trigger_to_use = "HLT_Jet70U";
   }
   else if (hardest_pt_value > 114) {
      trigger_to_use = "HLT_Jet50U";
   }
   else if (hardest_pt_value > 84) {
      trigger_to_use = "HLT_Jet30U";
   }
   else if (hardest_pt_value > 56) {
      trigger_to_use = "HLT_Jet15U";
   }
   else if (hardest_pt_value > 37) {
      trigger_to_use = "HLT_L1Jet6U";
   }
   else {
      trigger_to_use = "HLT_MinBiasPixel_SingleTrack";
   }
   
   // Here, we just return the trigger that was supposed to fire, not caring whether it actually did or not.
   // A check on whether it actually fired or not will be done in the analysis.cc file itself.

   return trigger_to_use;
}

bool MODEvent::read_event(ifstream & data_file) {
   // event_being_read = MODEvent();

   string line;
   while(getline(data_file, line)) {
      istringstream iss(line);
      
      int event_number, run_number;
      string tag, run_keyword, event_keyword;

      iss >> tag;
      istringstream stream(line);
      if (tag == "BeginEvent") {
         stream >> tag >> run_keyword >> run_number >> event_keyword >> event_number;
         set_event_number(event_number);
         set_run_number(run_number);
         set_particles_trigger_type("PFC");
      }
      else if (tag == "PFC") {
         try {
            add_particle(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format!");
            cout << "Something went wrong!" << endl;
         }
      }
      else if (tag == "trig") {
         try {
            add_trigger(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format!");
            cout << "Something went wrong!" << endl;
         }
      }
      else if (tag == "EndEvent") {
         return true;
      }
   }

   return false;
}

ostream& operator<< (ostream& os, const MODEvent& event) {
   os << event.make_string();
   return os;
}
