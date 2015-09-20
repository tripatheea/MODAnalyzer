#include "event.h"

using namespace std;
using namespace fastjet;

MOD::Event::Event(int run_number, int event_number, int lumi_block, double inst_lumi) : _run_number(run_number), _event_number(event_number) {}

MOD::Event::Event(int run_number, int event_number, int version, std::pair<std::string, std::string> data_type, MOD::Condition condition, vector<MOD::Trigger> triggers, vector<MOD::PFCandidate> particles, vector<fastjet::PseudoJet> pseudojets, vector<MOD::CalibratedJet> CMS_jets, vector<fastjet::PseudoJet> CMS_pseudojets) : 
_run_number(run_number), _event_number(event_number),  _version(version), _data_type(data_type), _condition(condition), _triggers(triggers), _particles(particles), _pseudojets(pseudojets), _CMS_jets(CMS_jets), _CMS_pseudojets(CMS_pseudojets) 
{}



MOD::Event::Event() {}

int MOD::Event::event_number() const {
   return _condition.event_number();
}

int MOD::Event::run_number() const {
   return _condition.run_number();
}

int MOD::Event::version() const {
   return _version;
}

pair<string, string> MOD::Event::data_type() const {
   return _data_type;
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

// const vector<PseudoJet> MOD::Event::charged_pseudojets(double pt_cut) const {
   
//    vector<PseudoJet> charged_pseudojets;

//    for (unsigned i = 0; i < charged_particles().size(); i++) {
//       if (charged_particles()[i].pseudojet().pt() >= pt_cut) {
//          charged_pseudojets.push_back(charged_particles()[i].pseudojet());
//       }
//    }

//    return charged_pseudojets;

// }


const vector<MOD::PFCandidate> & MOD::Event::particles() const {
   return _particles;
}

// const vector<MOD::PFCandidate> MOD::Event::charged_particles() const {
   
//    vector<MOD::PFCandidate> charged_particles;

//    for (unsigned i = 0; i < _particles.size(); i++) {
//       if ( (abs(_particles[i].pdgId()) == 211) || (abs(_particles[i].pdgId()) == 11) || (abs(_particles[i].pdgId()) == 13) ) {
//          charged_particles.push_back(_particles[i]);
//       }
//    }
//    return charged_particles;
// }




const vector<MOD::CalibratedJet> & MOD::Event::CMS_jets() const {
   return _CMS_jets;
}

const vector<fastjet::PseudoJet> & MOD::Event::CMS_pseudojets() const {
   return _CMS_pseudojets;
}

const vector<fastjet::PseudoJet> & MOD::Event::fastjet_pseudojets() const {
   return _fastjet_pseudojets;
}


void MOD::Event::add_particle(istringstream & input_stream) {
   MOD::PFCandidate new_particle(input_stream);

   _particles.push_back(new_particle);
   _pseudojets.push_back(PseudoJet(new_particle.pseudojet()));
}

void MOD::Event::add_condition(istringstream & input_stream) {
   MOD::Condition new_condition(input_stream);

   _condition = new_condition;
}

void MOD::Event::add_CMS_jet(istringstream & input_stream) {
   MOD::CalibratedJet new_jet(input_stream);

   _CMS_jets.push_back(new_jet);
   _CMS_pseudojets.push_back(new_jet.uncorrected_pseudojet());      
   
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

const bool MOD::Event::trigger_exists(string trigger_name) const {
   return trigger_by_name(trigger_name).is_valid();
}

const vector<MOD::Trigger> & MOD::Event::triggers() const {
   return _triggers;
}

string MOD::Event::make_string() const {
   stringstream file_to_write;
   
   file_to_write << "BeginEvent Version " << _version << " " << _data_type.first << " " << _data_type.second << endl;
   

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
   for(unsigned int i = 0; i < _CMS_jets.size(); i++) {
      if (i == 0)
         file_to_write << _CMS_jets[i].make_header_string();
      file_to_write << _CMS_jets[i];
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


bool MOD::Event::read_event(istream & data_stream) {

   string line;
   while(getline(data_stream, line)) {
      istringstream iss(line);

      int version;
      string tag, version_keyword, a, b;

      iss >> tag;      
      istringstream stream(line);

      if (tag == "BeginEvent") {

         stream >> tag >> version_keyword >> version >> a >> b;
         
         // cout << "BeginEvent" << endl;

         set_version(version);
         set_data_type(a, b);
      }
      else if (tag == "PFC") {
         try {
            // cout << "PFC" << endl;
            add_particle(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format PFC! Something's wrong with the way PFCs have been written.");
         }
      }
      else if (tag == "AK5") {
         try {
            // cout << "AK5" << endl;
            add_CMS_jet(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format AK5! Something's wrong with the way jets have been written.");
         }
      }
      else if (tag == "Trig") {
         try {
            // cout << "Trig" << endl;
            add_trigger(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format! Something's wrong with the way triggers have been written.");
         }
      }
      else if (tag == "Cond") {
         try {
            // cout << "Cond" << endl;
            add_condition(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format! Something's wrong with the way conditions have been written.");
         }
      }
      else if (tag == "EndEvent") {
         // cout << "EndEvent" << endl;
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
   bool fired = _assigned_trigger.is_valid() && _assigned_trigger.fired();
   return fired;
}

int MOD::Event::assigned_trigger_prescale() const {
   return _assigned_trigger.prescale();
}

void MOD::Event::set_assigned_trigger() {
   
   string trigger_to_use = "";
   string trigger = "";

   // First, figure out the hardest pT.
   MOD::CalibratedJet trigger_jet = _trigger_jet;

   
   double hardest_pT = trigger_jet.corrected_pseudojet().pt();


   // Next, lookup which trigger to use based on the pT value of the hardest jet.

   if ( (hardest_pT > 325) && ( trigger_exists("HLT_Jet180U") || trigger_exists("HLT_Jet180U_v1") || trigger_exists("HLT_Jet180U_v2") || trigger_exists("HLT_Jet180U_v3") ) ) {
      trigger = "HLT_Jet180U";
   }
   else if ( (hardest_pT > 260) && ( trigger_exists("HLT_Jet140U") || trigger_exists("HLT_Jet140U_v1") || trigger_exists("HLT_Jet140U_v2") || trigger_exists("HLT_Jet140U_v3") ) ) {
      trigger = "HLT_Jet140U";
   }
   else if ( (hardest_pT > 196) && ( trigger_exists("HLT_Jet100U") || trigger_exists("HLT_Jet100U_v1") || trigger_exists("HLT_Jet100U_v2") || trigger_exists("HLT_Jet100U_v3") ) ) {
      trigger = "HLT_Jet100U";
   }
   else if ( (hardest_pT > 153) && ( trigger_exists("HLT_Jet70U") || trigger_exists("HLT_Jet70U_v1") || trigger_exists("HLT_Jet70U_v2") || trigger_exists("HLT_Jet70U_v3") ) ) {
      trigger = "HLT_Jet70U";
   }
   else if ( (hardest_pT > 114) && ( trigger_exists("HLT_Jet50U") || trigger_exists("HLT_Jet50U_v1") || trigger_exists("HLT_Jet50U_v2") || trigger_exists("HLT_Jet50U_v3") ) ) {
      trigger = "HLT_Jet50U";
   }
   else if ( (hardest_pT > 84) && ( trigger_exists("HLT_Jet30U") || trigger_exists("HLT_Jet30U_v1") || trigger_exists("HLT_Jet30U_v2") || trigger_exists("HLT_Jet30U_v3") ) ) {
      trigger = "HLT_Jet30U";
   }
   else if ( (hardest_pT > 56) && ( trigger_exists("HLT_Jet15U") || trigger_exists("HLT_Jet15U_v1") || trigger_exists("HLT_Jet15U_v2") || trigger_exists("HLT_Jet15U_v3") ) ) {
      trigger = "HLT_Jet15U";
   }
   else if ( (hardest_pT > 37) && ( trigger_exists("HLT_L1Jet6U") || trigger_exists("HLT_L1Jet6U_v1") || trigger_exists("HLT_L1Jet6U_v2") || trigger_exists("HLT_L1Jet6U_v3") ) ) {
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

   if (trigger_to_use != "") {
      _assigned_trigger_name = trigger_to_use;
      _assigned_trigger = trigger_by_name(trigger_to_use);   
   }
   else {
      _assigned_trigger_name = "";
      _assigned_trigger = Trigger();
   }
}


void MOD::Event::set_trigger_jet() {
   // Get hardest jet, apply JEC, and then eta cut.

   vector<MOD::CalibratedJet> CMS_jets = _CMS_jets;

   vector<MOD::CalibratedJet> processed_jets = apply_jet_energy_corrections(CMS_jets);
   processed_jets = apply_eta_cut(processed_jets, 2.4);


   // Then, sort the jets and store the hardest one as _trigger_jet.
   if (processed_jets.size() > 0) {
      sort(processed_jets.begin(), processed_jets.end());
      
      auto it = std::find(processed_jets.begin(), processed_jets.end(), processed_jets[0]);
      auto index = std::distance(processed_jets.begin(), it);
      _trigger_jet = _CMS_jets[index];
   }
   else {
      _trigger_jet = CalibratedJet();
   }
}


void MOD::Event::set_closest_fastjet_jet_to_trigger_jet() {

   if (_trigger_jet.is_valid()) {

      // Cluster PFCandidates using FastJet.
      JetDefinition jet_def(antikt_algorithm, 0.5);
      ClusterSequence cs(_pseudojets, jet_def);
      vector<PseudoJet> fastjet_jets = sorted_by_pt(cs.inclusive_jets(3.0));

      // Loop through all FastJet pseudojets, calculating delta R for each one.
      vector<double> delta_Rs;
      for (unsigned i = 0; i < fastjet_jets.size(); i++) {
         delta_Rs.push_back( _trigger_jet.uncorrected_pseudojet().delta_R(fastjet_jets[i]) );
      }
      
      // Find the index of the fastjet jet that has the lowest delta R. This will be the jet that's "closest" to the CMS jet. 
      int index = -1;
      double delta_R = std::numeric_limits<double>::max();
      for (unsigned i = 0; i < delta_Rs.size(); i++) {
         if (delta_Rs[i] <= delta_R) {
            delta_R = delta_Rs[i];
            index = i;
         }
      }

      if (index >= 0) {
         // We now have the corresponding "hardest" FastJet jet.
         _closest_fastjet_jet_to_trigger_jet = fastjet_jets[index];
         _closest_fastjet_jet_to_trigger_jet_constituents = fastjet_jets[index].constituents();
         return;   
      }
   }

   return;
}


const MOD::CalibratedJet MOD::Event::trigger_jet() const {
   return _trigger_jet;
}

void MOD::Event::set_trigger_jet_is_matched() {

   fastjet::PseudoJet trigger_fastjet = _trigger_jet.uncorrected_pseudojet();
   fastjet::PseudoJet closest_fastjet_jet_to_trigger_jet = _closest_fastjet_jet_to_trigger_jet;

   // HOW CAN YOU COMPARE THE NUMBER OF CONSTITUENTS WHEN CMS AK5 JETS HAVE NO CONSTITUENTS?

   // // Compare the number of constituents first.
   // if (trigger_fastjet.constituents().size() != closest_fastjet_jet_to_trigger_jet.constituents().size()) {
   //    _trigger_jet_is_matched = false;
   //    return;
   // }

   // Next, compare if the 4-vector matches upto 10e-4 precision or not.
   double tolerance = pow(10, -3);
   if ( ( abs(trigger_fastjet.px() - closest_fastjet_jet_to_trigger_jet.px()) < tolerance ) && ( abs(trigger_fastjet.py() - closest_fastjet_jet_to_trigger_jet.py()) < tolerance ) && ( abs(trigger_fastjet.pz() - closest_fastjet_jet_to_trigger_jet.pz()) < tolerance ) && ( abs(trigger_fastjet.E() - closest_fastjet_jet_to_trigger_jet.E()) < tolerance ) ) {      
      _trigger_jet_is_matched = true;
      return;
   }

   _trigger_jet_is_matched = false;
   return;
}


fastjet::PseudoJet MOD::Event::closest_fastjet_jet_to_trigger_jet() {
   return _closest_fastjet_jet_to_trigger_jet;
}


std::vector<fastjet::PseudoJet> MOD::Event::closest_fastjet_jet_to_trigger_jet_constituents() {
   return _closest_fastjet_jet_to_trigger_jet_constituents;
}

void MOD::Event::establish_properties() {
   
   // Cluster PFCandidates into AK5 Jets using FastJet.
   vector<PseudoJet> pfcandidates = pseudojets();

   JetDefinition jet_def(antikt_algorithm, 0.5);
   ClusterSequence cs(pseudojets(), jet_def);
   vector<PseudoJet> ak5_jets = sorted_by_pt(cs.inclusive_jets(3.0));
   _fastjet_pseudojets = ak5_jets;

   // cout << "trigger_jet" << endl;
   // First of all, assign _trigger_jet.
   set_trigger_jet();

   // cout << "closest fastjet" << endl;
   // Next, find out the specific FastJet that's closest to _trigger_jet.
   set_closest_fastjet_jet_to_trigger_jet();

   // cout << "trigger_jet is matched" << endl;
   set_trigger_jet_is_matched();

   // cout << "assigned trigger" << endl;
   set_assigned_trigger();   

}





// MOD::PFCandidate MOD::Event::hardest_pfcandidate() {
//    // Get PFCandidates.
//    vector<MOD::PFCandidate> particles = _particles;

//    if (particles.size() > 0) {
//       // Sort by pT.
//       sort(particles.begin(), particles.end());

//       // Return the first PFCandidate.
//       return particles[0];
//    }
//    else {
//       throw runtime_error("No PFCandidate found!");
//    }
// }



vector<MOD::CalibratedJet> MOD::Event::apply_jet_energy_corrections(vector<MOD::CalibratedJet> jets) const {

   vector<MOD::CalibratedJet> jec_corrected_jets;

   for (unsigned i = 0; i < jets.size(); i++) {
      jec_corrected_jets.push_back(jets[i].corrected_jet());
   }

   return jec_corrected_jets;
}

vector<MOD::CalibratedJet> MOD::Event::apply_eta_cut(vector<MOD::CalibratedJet> jets, double eta_cut) const {

   vector<MOD::CalibratedJet> jets_with_eta_cut;

   for (unsigned i = 0; i < jets.size(); i++) {
      if ( abs(jets[i].uncorrected_pseudojet().eta()) < eta_cut ) {
         jets_with_eta_cut.push_back(jets[i]);
      }
   }

   return jets_with_eta_cut;
}

bool MOD::Event::trigger_jet_is_matched() const {
   return _trigger_jet_is_matched;
}



namespace MOD {
   ostream& operator<< (ostream& os, const Event& event) {
      os << event.make_string();
      return os;
   }
}
