#include "event.h"

using namespace std;
using namespace fastjet;

MOD::Event::Event(int run_number, int Event_number, int lumi_block, double inst_lumi) : _run_number(run_number), _event_number(Event_number), _hardest_pt(std::numeric_limits<double>::max()), _lumi_block(lumi_block), _inst_lumi(inst_lumi) {}

MOD::Event::Event() :  _hardest_pt(std::numeric_limits<double>::max()) {}

int MOD::Event::event_number() const {
   return _event_number;
}

int MOD::Event::run_number() const {
   return _run_number;
}

int MOD::Event::lumi_block() const {
   return _lumi_block;
}

double MOD::Event::inst_lumi() const {
   return _inst_lumi;
}

void MOD::Event::set_run_number(int run_number) {
   _run_number = run_number;
}

void MOD::Event::set_event_number(int event_number) {
   _event_number = event_number;
}

void MOD::Event::set_lumi_block(int lumi_block) {
   _lumi_block = lumi_block;
}

void MOD::Event::set_inst_lumi(double inst_lumi) {
   _inst_lumi = inst_lumi;
}

const vector<PseudoJet> & MOD::Event::pseudojets() const {
   return _pseudojets;
}

const vector<PseudoJet> & MOD::Event::calibrated_pseudojets() const {
   return _calibrated_pseudojets;
}

const vector<MOD::CalibratedJet> & MOD::Event::calibrated_jets() const {
   return _calibrated_jets;
}



const vector<MOD::PFCandidate> & MOD::Event::particles() const {
   return _particles;
}

const vector<MOD::CalibratedJet> MOD::Event::corrected_calibrated_jets() const {
   
   vector<MOD::CalibratedJet> corrected_jets;
   vector<MOD::CalibratedJet> jets = calibrated_jets();

   for (unsigned i = 0; i < jets.size(); i++) {
      corrected_jets.push_back(jets[i].corrected_jet());
   }

   return corrected_jets;
}


const vector<fastjet::PseudoJet> MOD::Event::corrected_calibrated_pseudojets() const {
   vector<fastjet::PseudoJet> corrected_psuedojets;
   vector<MOD::CalibratedJet> jets = calibrated_jets();

   for (unsigned i = 0; i < jets.size(); i++) {
      corrected_psuedojets.push_back(jets[i].corrected_jet().pseudojet());
   }

   return corrected_psuedojets;
}


void MOD::Event::add_particle(istringstream & input_stream) {
   MOD::PFCandidate new_particle(input_stream);

   _particles.push_back(new_particle);
   _pseudojets.push_back(PseudoJet(new_particle.pseudojet()));
}

void MOD::Event::add_calibrated_jet(istringstream & input_stream) {
   MOD::CalibratedJet new_jet(input_stream);

   _calibrated_jets.push_back(new_jet);
   _calibrated_pseudojets.push_back(PseudoJet(new_jet.pseudojet()));      
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
   
   file_to_write << "BeginEvent Run " << _run_number << " Event " << _event_number << " LumiSection " << _lumi_block << " AvgInstLumi " << _inst_lumi << endl;
   
   // First, write out all triggers.
   
   for(unsigned int i = 0; i < _triggers.size(); i++) {
      if (i == 0)
         file_to_write << _triggers[i].make_header_string();
      file_to_write << _triggers[i];
   }

   // Next, write out AK5 calibrated jets.
   for(unsigned int i = 0; i < _calibrated_jets.size(); i++) {
      if (i == 0)
         file_to_write << _calibrated_jets[i].make_header_string();
      file_to_write << _calibrated_jets[i];
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

double MOD::Event::hardest_pt() const {
   return _hardest_pt;
}

bool MOD::Event::read_event(istream & data_stream) {
   string line;
   while(getline(data_stream, line)) {
      istringstream iss(line);
      
      int event_number, run_number, lumi_section, event_serial_number;
      string tag, run_keyword, event_keyword, lumi_section_keyword, inst_lumi_keyword, event_serial_number_keyword;
      double inst_lumi;

      iss >> tag;      
      istringstream stream(line);
      if (tag == "BeginEvent") {
         stream >> tag >> run_keyword >> run_number >> event_keyword >> event_number >> lumi_section_keyword >> lumi_section >> inst_lumi_keyword >> inst_lumi >> event_serial_number_keyword >> event_serial_number;
         set_event_number(event_number);
         set_run_number(run_number);
         set_lumi_block(lumi_section);
         set_inst_lumi(inst_lumi);
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
   double hardest_pt_value = hardest_pt();

   if (hardest_pt_value == std::numeric_limits<double>::max()) {
      throw runtime_error("You need to set _trigger_hardest_pt before trying to retrieve assigned_trigger_name.");
   }

   // Next, lookup which trigger to use based on the pt value of the hardest jet.

   string trigger_to_use;
   if (hardest_pt_value > 153) {
      trigger_to_use = "HLT_Jet70U_v3";
   }
   else if (hardest_pt_value > 114) {
      trigger_to_use = "HLT_Jet50U_v3";
   }
   else if (hardest_pt_value > 84) {
      trigger_to_use = "HLT_Jet30U_v3";
   }
   else if (hardest_pt_value > 56) {
      trigger_to_use = "HLT_Jet15U_v3";
   }
   else if (hardest_pt_value > 37) {
      trigger_to_use = "HLT_L1Jet6U_v3";
   }
   else {
      trigger_to_use = "HLT_MinBiasPixel_SingleTrack";
   }


   _assigned_trigger_name = trigger_to_use;
   _assigned_trigger = trigger_by_name(trigger_to_use);
}

void MOD::Event::set_hardest_pt() {
   // Set the hardest pt of the AK5 jets.

   // Just use the jets we read from the data file.

   vector<PseudoJet> clustered_jets = _calibrated_pseudojets;
   
   double hardest_pt_value = 0.0;
   for (unsigned int i = 0; i < clustered_jets.size(); i++) {
      if (hardest_pt_value < clustered_jets[i].pt()) {
         hardest_pt_value = clustered_jets[i].pt();
      }
   }

   _hardest_pt = hardest_pt_value;

}

void MOD::Event::establish_properties() {
   set_hardest_pt();
   set_assigned_trigger();
}

double MOD::Event::hardest_jet_JEC() {
   // Get CMS Jets.
   vector<MOD::CalibratedJet> cms_jets = _calibrated_jets;
   
   if (cms_jets.size() > 0) {
      // Sort by pt.
      sort(cms_jets.begin(), cms_jets.end());

      // Return JEC of the first element.
      return cms_jets[0].JEC();
   }

   return cms_jets[0].JEC();
}

namespace MOD {
   ostream& operator<< (ostream& os, const Event& event) {
      os << event.make_string();
      return os;
   }
}
