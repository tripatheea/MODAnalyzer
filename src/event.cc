#include "../interface/event.h"

using namespace std;
using namespace fastjet;

MOD::Event::Event(int run_number, int Event_number) : _run_number(run_number), _event_number(Event_number) {}

MOD::Event::Event() {}

int MOD::Event::event_number() const {
   return _event_number;
}

int MOD::Event::run_number() const {
   return _run_number;
}

void MOD::Event::set_run_number(int run_number) {
   _run_number = run_number;
}

void MOD::Event::set_event_number(int event_number) {
   _event_number = event_number;
}

const vector<PseudoJet> & MOD::Event::pseudojets() const {
   return _pseudojets;
}

const vector<PseudoJet> & MOD::Event::calibrated_jets_pseudojets_ak5() const {
   return _calibrated_jets_pseudojets_ak5;
}

const vector<MOD::CalibratedJet> & MOD::Event::calibrated_jets_ak5() const {
   return _calibrated_jets_ak5;
}

const vector<PseudoJet> & MOD::Event::calibrated_jets_pseudojets_ak7() const {
   return _calibrated_jets_pseudojets_ak7;
}

const vector<MOD::CalibratedJet> & MOD::Event::calibrated_jets_ak7() const {
   return _calibrated_jets_ak7;
}

const vector<MOD::PFCandidate> & MOD::Event::particles() const {
   return _particles;
}

void MOD::Event::add_particle(istringstream & input_stream) {
   MOD::PFCandidate new_particle = MOD::PFCandidate(input_stream);
   _particles.push_back(new_particle);
   _pseudojets.push_back(PseudoJet(new_particle.pseudojet()));
}

void MOD::Event::add_calibrated_jet(istringstream & input_stream) {
   MOD::CalibratedJet new_jet = MOD::CalibratedJet(input_stream);
   
   if (new_jet.algorithm() == "AK5") {
      _calibrated_jets_ak5.push_back(new_jet);
      _calibrated_jets_pseudojets_ak5.push_back(PseudoJet(new_jet.pseudojet()));   
   }
   else if (new_jet.algorithm() == "AK7") {
      _calibrated_jets_ak7.push_back(new_jet);
      _calibrated_jets_pseudojets_ak7.push_back(PseudoJet(new_jet.pseudojet()));   
   }
   
}

void MOD::Event::add_trigger(istringstream & input_stream) {
   _triggers.push_back(MOD::Trigger(input_stream));
}

const MOD::Trigger & MOD::Event::trigger_by_name(string name) const {
   for(int i = 0; i < triggers().size(); i++) {
      const MOD::Trigger& current_trigger = triggers().at(i);

      if (current_trigger.name() == name) {
         return current_trigger;
      }
   }

   MOD::Trigger * empty_trigger = new Trigger();
   return *empty_trigger;
}

const vector<MOD::Trigger> & MOD::Event::triggers() const {
   return _triggers;
}

string MOD::Event::make_string() const {
   stringstream file_to_write;
   
   file_to_write << "BeginEvent Run " << _run_number << " Event " << _event_number << endl;
   
   // First, write out all triggers.
   file_to_write << _triggers[0].make_header_string();
   for(int i = 0; i < _triggers.size(); i++) {
      file_to_write << _triggers[i];
   }

   // Next, write out AK5 calibrated jets.
   file_to_write << _calibrated_jets_ak5[0].make_header_string();
   for(int i = 0; i < _calibrated_jets_ak5.size(); i++) {
      file_to_write << _calibrated_jets_ak5[i];
   }

   // AK7 calibrated jets.
   file_to_write << _calibrated_jets_ak7[0].make_header_string();
   for(int i = 0; i < _calibrated_jets_ak7.size(); i++) {
      file_to_write << _calibrated_jets_ak7[i];
   }
   
   // Finally, write out all particles.
   file_to_write << _particles[0].make_header_string();
   for (int i = 0; i < _particles.size(); i++) {
      file_to_write << _particles[i];
   }

   file_to_write << "EndEvent" << endl;

   return file_to_write.str();
}

double MOD::Event::trigger_hardest_pt() const {
   return _trigger_hardest_pt;
}

string MOD::Event::assigned_trigger_name() const {

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

bool MOD::Event::read_event(ifstream & data_file) {
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
      }
      else if (tag == "PFC") {
         try {
            add_particle(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format PFC! Something's wrong with the way PFCs have been written.");
         }
      }
      else if ( (tag == "AK5") || (tag == "AK7") ) {
         try {
            add_calibrated_jet(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format AK! Something's wrong with the way jets have been written.");
         }
      }
      else if (tag == "trig") {
         try {
            add_trigger(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format! Something's wrong with the way triggers have been written.");
         }
      }
      else if (tag == "EndEvent") {
         establish_properties();
         return true;
      }
      else if (tag == "#") {
         // This line in the data file represents a comment. Just ignore it.
      }
      else {
         throw runtime_error("Invalid file format! Unrecognized tag!");
      }
   }

   return false;
}

const MOD::Trigger & MOD::Event::assigned_trigger() const {
   return _assigned_trigger;
}

bool MOD::Event::assigned_trigger_fired() const {
   return _assigned_trigger.fired();
}

int MOD::Event::assigned_trigger_prescale() const {
   return _assigned_trigger.prescale();
}

void MOD::Event::set_assigned_trigger() {
   _assigned_trigger = trigger_by_name(assigned_trigger_name());
}

void MOD::Event::set_trigger_hardest_pt() {
   // Get the hardest pt of the AK5 jets.

   // If we already have AK5 jets, use that, else run clustering first using FastJet.

   vector<PseudoJet> clustered_jets;
   if (_calibrated_jets_ak5.size() != 0) {
      // Just use the jets we read from the data file.
      clustered_jets = _calibrated_jets_pseudojets_ak5;
   }
   else {
      // Run the clustering, extract the jets using fastjet.
      JetDefinition jet_def(antikt_algorithm, 0.5);
      ClusterSequence cs(pseudojets(), jet_def);
      clustered_jets = cs.inclusive_jets(0.0);   
   }

   double hardest_pt = 0.0;
   for (unsigned int i = 0; i < clustered_jets.size(); i++) {
      if (hardest_pt < clustered_jets[i].pt()) {
         hardest_pt = clustered_jets[i].pt();
      }
   }

   _trigger_hardest_pt = hardest_pt;
}

void MOD::Event::establish_properties() {
   set_assigned_trigger();
   set_trigger_hardest_pt();
}


namespace MOD {
   ostream& operator<< (ostream& os, const Event& event) {
      os << event.make_string();
      return os;
   }
}
