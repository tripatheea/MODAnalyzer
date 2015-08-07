#include "event.h"

using namespace std;
using namespace fastjet;

MOD::Event::Event(int run_number, int Event_number, int lumi_block, double inst_lumi) : _run_number(run_number), _event_number(Event_number), _hardest_corrected_pt(std::numeric_limits<double>::max()), _hardest_uncorrected_pt(std::numeric_limits<double>::max()) {}

MOD::Event::Event() :  _hardest_corrected_pt(std::numeric_limits<double>::max()), _hardest_uncorrected_pt(std::numeric_limits<double>::max()) {}

int MOD::Event::event_number() const {
   return _event_number;
}

int MOD::Event::run_number() const {
   return _run_number;
}

int MOD::Event::version() const {
   return _version;
}

pair<string, string> MOD::Event::data_type() const {
   return _data_type;
}

void MOD::Event::set_run_number(int run_number) {
   _run_number = run_number;
}

void MOD::Event::set_event_number(int event_number) {
   _event_number = event_number;
}

void MOD::Event::set_version(int version) {
   _version = version;
}

void MOD::Event::set_data_type(string a, string b) {
   _data_type = make_pair(a, b);;
}

const vector<PseudoJet> MOD::Event::pseudojets(double pt_cut) const {
   
   vector<PseudoJet> pfcandidates;

   for (unsigned i = 0; i < _pseudojets.size(); i++) {
      if (_pseudojets[i].pt() >= pt_cut) {
         pfcandidates.push_back(_pseudojets[i]);
      }
   }

   return pfcandidates;

}

const vector<PseudoJet> MOD::Event::charged_pseudojets(double pt_cut) const {
   
   vector<PseudoJet> charged_pseudojets;

   for (unsigned i = 0; i < charged_particles().size(); i++) {
      if (charged_particles()[i].pseudojet().pt() >= pt_cut) {
         charged_pseudojets.push_back(charged_particles()[i].pseudojet());
      }
   }

   return charged_pseudojets;

}


const vector<MOD::PFCandidate> & MOD::Event::particles() const {
   return _particles;
}

const vector<MOD::PFCandidate> MOD::Event::charged_particles() const {
   
   vector<MOD::PFCandidate> charged_particles;

   for (unsigned i = 0; i < _particles.size(); i++) {
      if ( (abs(_particles[i].pdgId()) == 211) || (abs(_particles[i].pdgId()) == 11) || (abs(_particles[i].pdgId()) == 13) ) {
         charged_particles.push_back(_particles[i]);
      }
   }
   return charged_particles;
}

const vector<MOD::CalibratedJet> MOD::Event::corrected_calibrated_jets() const {
   
   vector<MOD::CalibratedJet> corrected_jets;
   vector<MOD::CalibratedJet> jets = corrected_calibrated_jets();

   for (unsigned i = 0; i < jets.size(); i++) {
      corrected_jets.push_back(jets[i].corrected_jet());
   }

   return corrected_jets;
}

const vector<fastjet::PseudoJet> MOD::Event::corrected_calibrated_pseudojets() const {
   vector<fastjet::PseudoJet> corrected_psuedojets;
   vector<MOD::CalibratedJet> jets = corrected_calibrated_jets();

   for (unsigned i = 0; i < jets.size(); i++) {
      corrected_psuedojets.push_back(jets[i].corrected_jet().pseudojet());
   }

   return corrected_psuedojets;
}

const vector<PseudoJet> & MOD::Event::uncorrected_calibrated_pseudojets() const {
   return _uncorrected_calibrated_pseudojets;
}

const vector<MOD::CalibratedJet> & MOD::Event::uncorrected_calibrated_jets() const {
   return _uncorrected_calibrated_jets;
}


void MOD::Event::add_particle(istringstream & input_stream) {
   MOD::PFCandidate new_particle(input_stream);

   _particles.push_back(new_particle);
   _pseudojets.push_back(PseudoJet(new_particle.pseudojet()));
}

void MOD::Event::add_conditions(istringstream & input_stream) {
   MOD::Condition new_condition(input_stream);

   _condition = new_condition;
}

void MOD::Event::add_calibrated_jet(istringstream & input_stream) {
   MOD::CalibratedJet new_jet(input_stream);

   _uncorrected_calibrated_jets.push_back(new_jet);
   _uncorrected_calibrated_pseudojets.push_back(new_jet.pseudojet());      
   
   _corrected_calibrated_jets.push_back(new_jet.corrected_jet());
   _corrected_calibrated_pseudojets.push_back(new_jet.JEC() * new_jet.pseudojet());      
}

void MOD::Event::add_trigger(istringstream & input_stream) {
   _triggers.push_back(MOD::Trigger(input_stream));
}

const MOD::Trigger MOD::Event::trigger_by_name(string name) const {
   for(unsigned int i = 0; i < triggers().size(); i++) {
      const MOD::Trigger& current_trigger = triggers().at(i);

      if (current_trigger.name() == name) {
         return current_trigger;
      }
   }

   return Trigger();
}


const vector<MOD::Trigger> & MOD::Event::triggers() const {
   return _triggers;
}

string MOD::Event::make_string() const {
   stringstream file_to_write;
   
   file_to_write << "BeginEvent Version " << _version << " " << _data_type.first << " " << _data_type.second << " Run " << _run_number << " Event " << _event_number << endl;
   

   // First, write out conditions.
   file_to_write << _condition.make_header_string();
   file_to_write << _condition;

   // Then, write out all triggers.
   
   for(unsigned int i = 0; i < _triggers.size(); i++) {
      if (i == 0)
         file_to_write << _triggers[i].make_header_string();
      file_to_write << _triggers[i];
   }

   // Next, write out AK5 calibrated jets.
   for(unsigned int i = 0; i < _uncorrected_calibrated_jets.size(); i++) {
      if (i == 0)
         file_to_write << _uncorrected_calibrated_jets[i].make_header_string();
      file_to_write << _uncorrected_calibrated_jets[i];
   }
   
   // Finally, write out all particles.
   for (unsigned int i = 0; i < _particles.size(); i++) {
      if (i == 0)
         file_to_write << _particles[i].make_header_string();
      file_to_write << _particles[i];
   }

   file_to_write << "EndEvent" << endl;

   return file_to_write.str();
}

double MOD::Event::hardest_corrected_pt() const {
   return _hardest_corrected_pt;
}

double MOD::Event::hardest_uncorrected_pt() const {
   return _hardest_uncorrected_pt;
}

bool MOD::Event::read_event(istream & data_stream) {
   string line;
   while(getline(data_stream, line)) {
      istringstream iss(line);

      int event_number, run_number, version;
      string tag, run_keyword, event_keyword, version_keyword, a, b;

      iss >> tag;      
      istringstream stream(line);

      if (tag == "BeginEvent") {

         stream >> tag >> version_keyword >> version >> a >> b >> run_keyword >> run_number >> event_keyword >> event_number;
         
         set_event_number(event_number);
         set_run_number(run_number);
         set_version(version);
         set_data_type(a, b);
      }
      else if (tag == "PFC") {
         try {
            add_particle(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format PFC! Something's wrong with the way PFCs have been written.");
         }
      }
      else if (tag == "AK5") {
         try {
            add_calibrated_jet(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format AK5! Something's wrong with the way jets have been written.");
         }
      }
      else if (tag == "Trig") {
         try {
            add_trigger(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format! Something's wrong with the way triggers have been written.");
         }
      }
      else if (tag == "Cond") {
         try {
            add_conditions(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format! Something's wrong with the way conditions have been written.");
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
         throw runtime_error("Invalid file format! Unrecognized tag: " + tag + "!");
      }
   }

   return false;
}

string MOD::Event::assigned_trigger_name() const {
   return _assigned_trigger_name;
}

const MOD::Trigger MOD::Event::assigned_trigger() const {
   return _assigned_trigger;
}

bool MOD::Event::assigned_trigger_fired() const {
   return _assigned_trigger.fired();
}

int MOD::Event::assigned_trigger_prescale() const {
   return _assigned_trigger.prescale();
}

void MOD::Event::set_assigned_trigger() {
   double hardest_pt_value = hardest_corrected_pt();

   if (hardest_pt_value == std::numeric_limits<double>::max()) {
      throw runtime_error("You need to set _trigger_hardest_pt before trying to retrieve assigned_trigger_name.");
   }

   // Next, lookup which trigger to use based on the pt value of the hardest jet.

   
   string trigger_to_use;
   string trigger;

   if (hardest_pt_value > 153) {
      trigger = "HLT_Jet70U";
   }
   else if (hardest_pt_value > 114) {
      trigger = "HLT_Jet50U";
   }
   else if (hardest_pt_value > 84) {
      trigger = "HLT_Jet30U";
   }
   else if (hardest_pt_value > 56) {
      trigger = "HLT_Jet15U";
   }
   else if (hardest_pt_value > 37) {
      trigger = "HLT_L1Jet6U";
   }

   // Since there are multiple trigger versions, keep trying until you find the right one.

   if ( trigger_by_name(trigger).is_valid()) {
      trigger_to_use = trigger;
   }
   else if (trigger_by_name(trigger + "_v1").is_valid()) {
      trigger_to_use = trigger + "_v1";
   }
   else if (trigger_by_name(trigger + "_v2").is_valid()) {
      trigger_to_use = trigger + "_v2";
   }
   else if (trigger_by_name(trigger + "_v3").is_valid()) {
      trigger_to_use = trigger + "_v3";
   }



   _assigned_trigger_name = trigger_to_use;
   _assigned_trigger = trigger_by_name(trigger_to_use);
}

void MOD::Event::set_hardest_pt() {
   // Set the hardest pt of the AK5 jets.

   // Just use the jets we read from the data file.

   vector<PseudoJet> corrected_clustered_jets = _corrected_calibrated_pseudojets;
   
   double corrected_hardest_pt_value = 0.0;
   for (unsigned int i = 0; i < corrected_clustered_jets.size(); i++) {
      if (corrected_hardest_pt_value < corrected_clustered_jets[i].pt()) {
         corrected_hardest_pt_value = corrected_clustered_jets[i].pt();
      }
   }

   _hardest_corrected_pt = corrected_hardest_pt_value;

   // Uncorrected Jets.

   vector<PseudoJet> uncorrected_clustered_jets = _uncorrected_calibrated_pseudojets;
   
   double uncorrected_hardest_pt_value = 0.0;
   for (unsigned int i = 0; i < uncorrected_clustered_jets.size(); i++) {
      if (uncorrected_hardest_pt_value < uncorrected_clustered_jets[i].pt()) {
         uncorrected_hardest_pt_value = uncorrected_clustered_jets[i].pt();
      }
   }

   _hardest_uncorrected_pt = uncorrected_hardest_pt_value;

}

void MOD::Event::establish_properties() {
   set_hardest_pt();
   set_assigned_trigger();
}



MOD::CalibratedJet MOD::Event::hardest_corrected_jet() {
   // Get CMS Jets.
   vector<MOD::CalibratedJet> cms_jets = _corrected_calibrated_jets;
   
   if (cms_jets.size() > 0) {
      // Sort by pT.
      sort(cms_jets.begin(), cms_jets.end());

      // Return the first element.
      return cms_jets[0];
   }
   else {
      // throw runtime_error("No jet found!");
      return MOD::CalibratedJet();
   }
}

MOD::CalibratedJet MOD::Event::hardest_uncorrected_jet() {
   // Get CMS Jets.
   vector<MOD::CalibratedJet> cms_jets = _uncorrected_calibrated_jets;
   
   if (cms_jets.size() > 0) {
      // Sort by pT.
      sort(cms_jets.begin(), cms_jets.end());

      // Return the first element.
      return cms_jets[0];
   }
   else {
      // throw runtime_error("No jet found!");
      return MOD::CalibratedJet();
   }
}

MOD::PFCandidate MOD::Event::hardest_pfcandidate() {
   // Get PFCandidates.
   vector<MOD::PFCandidate> particles = _particles;

   if (particles.size() > 0) {
      // Sort by pT.
      sort(particles.begin(), particles.end());

      // Return the first PFCandidate.
      return particles[0];
   }
   else {
      throw runtime_error("No PFCandidate found!");
   }
}

namespace MOD {
   ostream& operator<< (ostream& os, const Event& event) {
      os << event.make_string();
      return os;
   }
}
