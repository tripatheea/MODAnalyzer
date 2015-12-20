#include "event.h"


using namespace std;
using namespace fastjet;

MOD::Event::Event(int run_number, int event_number, int lumi_block, double inst_lumi) : _run_number(run_number), _event_number(event_number) {}

// MOD::Event::Event(int run_number, int event_number, int version, std::pair<std::string, std::string> data_type, MOD::Condition condition, vector<MOD::Trigger> triggers, vector<MOD::PFCandidate> particles, vector<MOD::CalibratedJet> CMS_jets) : 
// _run_number(run_number), _event_number(event_number),  _version(version), _data_type(data_type), _condition(condition), _triggers(triggers), _particles(particles), _CMS_jets(CMS_jets) 
// {}



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

const vector<MOD::PFCandidate> & MOD::Event::particles() const {
   return _particles;
}


const vector<MOD::CalibratedJet> & MOD::Event::CMS_jets() const {
   return _CMS_jets;
}


const vector<fastjet::PseudoJet> & MOD::Event::fastjet_clustered_pseudojets() const {
   return _fastjet_clustered_pseudojets;
}


void MOD::Event::add_particle(istringstream & input_stream) {
   MOD::PFCandidate new_particle(input_stream);
   _particles.push_back(new_particle);
}

void MOD::Event::add_mc_truth_particle(istringstream & input_stream) {
   MOD::MCPFCandidate new_particle(input_stream);
   _mc_truth_particles.push_back(new_particle);
}

void MOD::Event::add_mc_reco_particle(istringstream & input_stream) {
   MOD::MCPFCandidate new_particle(input_stream);
   _mc_reco_particles.push_back(new_particle);
}

void MOD::Event::add_condition(istringstream & input_stream) {
   MOD::Condition new_condition(input_stream);
   _condition = new_condition;
}

void MOD::Event::add_CMS_jet(istringstream & input_stream) {
   MOD::CalibratedJet new_jet(input_stream);
   _CMS_jets.push_back(new_jet);
}



void MOD::Event::add_trigger(istringstream & input_stream) {
   _triggers.push_back(MOD::Trigger(input_stream));
}

const int MOD::Event::data_source() const {
   return _data_source;
}

void MOD::Event::set_data_source(int data_source) {
   // 0 => Experiment, 1 => MC_TRUTH, 2 => MC_RECO.
   _data_source = static_cast<data_source_t>(data_source);
}

const MOD::Trigger MOD::Event::trigger_by_name(string name) const {
   for (unsigned i = 0; i < triggers().size(); i++) {
      const MOD::Trigger& current_trigger = triggers().at(i);

      if (current_trigger.name() == name) 
         return current_trigger;
   }

   return Trigger();
}

const MOD::Trigger MOD::Event::trigger_by_short_name(string short_name) const {
   for (unsigned i = 0; i < triggers().size(); i++) {
      const MOD::Trigger& current_trigger = triggers().at(i);

      if (current_trigger.short_name() == short_name)
         return current_trigger;
   }

   return Trigger();
}

// You can give a trigger's full name or short name here.
const bool MOD::Event::trigger_exists(string trigger_name) const {
   return ( trigger_by_name(trigger_name).is_valid() || trigger_by_short_name(trigger_name).is_valid() );
}

const vector<MOD::Trigger> & MOD::Event::triggers() const {
   return _triggers;
}



const int MOD::Event::prescale() const {
   return _prescale;
}

void MOD::Event::set_prescale(int prescale) {
   _prescale = prescale;
}

string MOD::Event::make_string() const {
   stringstream file_to_write;
   
   int data_source = _data_source;  // EXPERIMENT = 0, MC_TRUTH = 1, MC_RECO = 2 

   
   file_to_write << "BeginEvent Version " << _version << " " << _data_type.first << " " << _data_type.second << " Prescale " << _prescale << endl;      
  

   if (data_source == 0) { // Data is from experiment.
      
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

   }
   else if (data_source == 1) {
      
      for (unsigned i = 0; i < _mc_truth_particles.size(); i++) {
         if (i == 0)
            file_to_write << _mc_truth_particles[i].make_header_string();
         file_to_write << _mc_truth_particles[i];
      }

   }
   else if (data_source == 2) {

      for (unsigned i = 0; i < _mc_reco_particles.size(); i++) {
         if (i == 0)
            file_to_write << _mc_reco_particles[i].make_header_string();
         file_to_write << _mc_reco_particles[i];
      }

   }
   else if (data_source == 3) {

      // This output is for pristine form. Data in "pristine form" is used as-it-is i.e. without any consideration of things like JEC.
      // What this means here is that we need to filter / apply all corrections here itself before outputing the event.

      // We only output the constituents of the hardest jet.
      // The reason for that is, we want to do our analysis on FastJet-clustered jets instead of on CMS jets. However, we need to apply JEC to our AK5 jets. 
      // Those JEC factors are known for CMS jets but not for FastJet-clustered jets. This means, we need to figure out a correspondance between CMS and FastJet-clustered jets.
      // This is hard to do for jets other than the hardest one.
      // So we select the hardest CMS jet, use delta R to find the corresponding FastJet-clustered jet and then just output constituents of that jet.

      // While it's possible to do that for all jets, we output just the hardest jet's constituents because all analyses are going to be performed on the hardest jet anyway.

      // PseudoJet closest_fastjet_jet_to_trigger_jet = closest_fastjet_jet_to_trigger_jet();

      file_to_write << "# PDPFC" << "              px              py              pz          energy   pdgId" << endl;

      cout << _pristine_particles.size() << endl;

      for (unsigned i = 0; i < _pristine_particles.size(); i++) {
         
         file_to_write << "  PDPFC"
                       << setw(16) << fixed << setprecision(8) << _pristine_particles[i].pseudojet().px()
                       << setw(16) << fixed << setprecision(8) << _pristine_particles[i].pseudojet().py()
                       << setw(16) << fixed << setprecision(8) << _pristine_particles[i].pseudojet().pz()
                       << setw(16) << fixed << setprecision(8) << _pristine_particles[i].pseudojet().E()
                       << setw(8) << noshowpos << _pristine_particles[i].pseudojet().user_index()
                       << endl;
      }

   }
   else {
      throw runtime_error("Invalid data source!");
   }

   
   file_to_write << "EndEvent" << endl;

   return file_to_write.str();
}


const MOD::MCCalibratedJet MOD::Event::hardest_mc_truth_jet() const {
   return _hardest_mc_truth_jet;
}

const MOD::MCCalibratedJet MOD::Event::hardest_mc_reco_jet() const {
   return _hardest_mc_reco_jet;
}


const MOD::PDCalibratedJet MOD::Event::hardest_pristine_jet() const {
   return _pristine_jets[0];
}

const fastjet::PseudoJet MOD::Event::hardest_jet() const {
   
   // EXPERIMENT = 0, MC_TRUTH = 1, MC_RECO = 2, PRISTINE = 3 
   
   if (_data_source == 0) {
      return _closest_fastjet_jet_to_trigger_jet;
   }
   else if (_data_source == 1) {
      return _hardest_mc_truth_jet.pseudojet();
   }
   else if (_data_source == 2) {
      return _hardest_mc_reco_jet.pseudojet();
   }
   else if (_data_source == 3) {
      return _pristine_jets[0].pseudojet();
   }

   return PseudoJet();

}


void MOD::Event::convert_to_pristine() {

   // Set the data source to "Pristine".
   _data_source = static_cast<data_source_t>(3);

   _prescale = _assigned_trigger.prescale();
   

   PseudoJet jec_corrected_jet = _closest_fastjet_jet_to_trigger_jet * _trigger_jet.JEC();
   vector<PseudoJet> jec_corrected_jet_constituents = jec_corrected_jet.constituents();
   MOD::PDCalibratedJet jet = PDCalibratedJet(jec_corrected_jet, "AK5");

   vector<MOD::PDPFCandidate> particles;

   for (unsigned i = 0; i < jec_corrected_jet_constituents.size(); i++) {
      particles.push_back( PDPFCandidate(jec_corrected_jet_constituents[i]) );
   }

   vector<MOD::PDCalibratedJet> jets{jet};

   _pristine_particles = particles;
   _pristine_jets = jets;

   // Empty pfcandidates, jets, triggers.
   _particles.clear();
   _CMS_jets.clear();
   _triggers.clear();
   _fastjet_clustered_pseudojets.clear();

}


void MOD::Event::set_hardest_truth_jet() {
   
   if (_mc_truth_jets.size() == 0) {
      _hardest_mc_truth_jet = MCCalibratedJet();
   }
   else {

      sort(_mc_truth_jets.begin(), _mc_truth_jets.end());

      _hardest_mc_truth_jet = _mc_truth_jets[0];

      // Recluster stuff to get the constituents.

      vector<PseudoJet> truth_particles_pseudojets = MOD::convert_to_pseudojets(_mc_truth_particles);

      JetDefinition jet_def(antikt_algorithm, 0.5);
      ClusterSequence cs(truth_particles_pseudojets, jet_def);
      vector<PseudoJet> fastjet_jets = sorted_by_pt(cs.inclusive_jets(0.0));
   }
}

void MOD::Event::set_hardest_reco_jet() {

   if (_mc_reco_jets.size() == 0) {
      _hardest_mc_reco_jet = MCCalibratedJet();
   }
   else {

      sort(_mc_reco_jets.begin(), _mc_reco_jets.end());
      _hardest_mc_reco_jet = _mc_reco_jets[0];

      // Recluster stuff to get the constituents.

      vector<PseudoJet> reco_particles_pseudojets = MOD::convert_to_pseudojets(_mc_reco_particles);

      JetDefinition jet_def(antikt_algorithm, 0.5);
      ClusterSequence cs(reco_particles_pseudojets, jet_def);
      vector<PseudoJet> fastjet_jets = sorted_by_pt(cs.inclusive_jets(0.0));
   }
}



bool MOD::Event::read_event(istream & data_stream) {

   string line;
   while(getline(data_stream, line)) {
      istringstream iss(line);

      int version, prescale;
      string tag, version_keyword, a, b, prescale_keyword;

      iss >> tag;      
      istringstream stream(line);

      if (tag == "BeginEvent") {

         stream >> tag >> version_keyword >> version >> a >> b >> prescale_keyword >> prescale;
         
         // cout << "BeginEvent" << endl;

         set_version(version);
         set_data_type(a, b);

         if ( (prescale_keyword != NULL) && (prescale != NULL)) {
            set_prescale(prescale);
         }
         else {
            set_prescale(0);
         }

      }
      else if (tag == "PFC") {
         try {
            // cout << "PFC" << endl;
            set_data_source(0);
            add_particle(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format PFC! Something's wrong with the way PFCs have been written.");
         }
      }
      else if (tag == "AK5") {
         try {
            // cout << "AK5" << endl;
            set_data_source(0);
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
      else if (tag == "TRUTH") {
         try {
            set_data_source(1);
            add_mc_truth_particle(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format! Something's wrong with the way Truth has been written. ;)");
         }
      }
      else if (tag == "RPFC") {
         try {
            set_data_source(2);
            add_mc_reco_particle(stream);
         }
         catch (exception& e) {
            throw runtime_error("Invalid file format! Something's wrong with the way RPFC has been written. ;)");
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
   // cout << _assigned_trigger.is_valid() << ", " << _assigned_trigger.fired() << endl;
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

   if ( (hardest_pT > 325) && trigger_exists("HLT_Jet180U") ) {
      trigger = "HLT_Jet180U";
   }
   else if ( (hardest_pT > 260) && trigger_exists("HLT_Jet140U") ) {
      trigger = "HLT_Jet140U";
   }
   else if ( (hardest_pT > 196) && trigger_exists("HLT_Jet100U") ) {
      trigger = "HLT_Jet100U";
   }
   else if ( (hardest_pT > 153) && trigger_exists("HLT_Jet70U") ) {
      trigger = "HLT_Jet70U";
   }
   else if ( (hardest_pT > 114) && trigger_exists("HLT_Jet50U") ) {
      trigger = "HLT_Jet50U";
   }
   else if ( (hardest_pT > 84) && trigger_exists("HLT_Jet30U") ) {
      trigger = "HLT_Jet30U";
   }
   else if ( (hardest_pT > 56) && trigger_exists("HLT_Jet15U") ) {
      trigger = "HLT_Jet15U";
   }
   else if ( (hardest_pT > 37) && trigger_exists("HLT_L1Jet6U") ) {
      trigger = "HLT_L1Jet6U";
   }

   trigger_to_use = trigger;

   if (trigger_to_use != "") {
      _assigned_trigger_name = trigger_to_use;
      _assigned_trigger = trigger_by_short_name(trigger_to_use);   
   }
   else {
      _assigned_trigger_name = "";
      _assigned_trigger = Trigger();
   }

}


void MOD::Event::set_trigger_jet() {
   // Get the hardest jet, apply JEC, and then eta cut.

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
      vector<PseudoJet> pseudojets = convert_to_pseudojets(_particles);


      JetDefinition jet_def(antikt_algorithm, 0.5);
      ClusterSequence * cs = new ClusterSequence(pseudojets, jet_def);
      vector<PseudoJet> fastjet_jets = sorted_by_pt(cs->inclusive_jets(3.0));

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

      // This preserves the scope of our ClusterSequence so that later internal structure (like constituents) of the jet can be retrieved.
      cs->delete_self_when_unused();

      if (index >= 0) {
         // We now have the corresponding "hardest" FastJet jet.
         _closest_fastjet_jet_to_trigger_jet = fastjet_jets[index];
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

   // Compare the number of constituents first.
   if ((unsigned) _trigger_jet.number_of_constituents() != (unsigned) _closest_fastjet_jet_to_trigger_jet.constituents().size()) {
      _trigger_jet_is_matched = false;
      return;
   }

   // Next, compare if the 4-vector matches upto 10e-4 precision or not.
   double tolerance = pow(10, -3);
   if ( ( abs(trigger_fastjet.px() - _closest_fastjet_jet_to_trigger_jet.px()) < tolerance ) && ( abs(trigger_fastjet.py() - _closest_fastjet_jet_to_trigger_jet.py()) < tolerance ) && ( abs(trigger_fastjet.pz() - _closest_fastjet_jet_to_trigger_jet.pz()) < tolerance ) && ( abs(trigger_fastjet.E() - _closest_fastjet_jet_to_trigger_jet.E()) < tolerance ) ) {      
      _trigger_jet_is_matched = true;
      return;
   }

   _trigger_jet_is_matched = false;
   return;
}


fastjet::PseudoJet MOD::Event::closest_fastjet_jet_to_trigger_jet() {
   return _closest_fastjet_jet_to_trigger_jet;
}



void MOD::Event::establish_properties() {

   JetDefinition jet_def(antikt_algorithm, 0.5);

   if (data_source() == 0) {  // Experiment 

      // Cluster PFCandidates into AK5 Jets using FastJet.
      vector<MOD::PFCandidate> pfcandidates = particles();
      vector<PseudoJet> pfcandidates_pseudojets = MOD::convert_to_pseudojets(pfcandidates);
   
      ClusterSequence cs(pfcandidates_pseudojets, jet_def);
      vector<PseudoJet> ak5_jets = sorted_by_pt(cs.inclusive_jets(3.0));
      _fastjet_clustered_pseudojets = ak5_jets;

      // First of all, assign _trigger_jet.
      set_trigger_jet();

      // Next, find out the specific FastJet that's closest to _trigger_jet.
      set_closest_fastjet_jet_to_trigger_jet();

      set_trigger_jet_is_matched();

      set_assigned_trigger();   
   }
   else if (data_source() == 1) {   // 1 => MC_TRUTH 
      
      // Recluster all truth particles to get "Truth Jets"
      
      vector<PseudoJet> truth_particles_pseudojets = convert_to_pseudojets(_mc_truth_particles);

      ClusterSequence cs(truth_particles_pseudojets, jet_def);
      vector<PseudoJet> truth_ak5_jets = sorted_by_pt(cs.inclusive_jets(0.0));

      // Create a vector of MOD::MCCalibratedJet.
      vector<MOD::MCCalibratedJet> truth_jets;
      for (unsigned i = 0; i < truth_ak5_jets.size(); i++) {
         truth_jets.push_back(MOD::MCCalibratedJet( truth_ak5_jets[i], "ak5" ));
      }

      _mc_truth_jets = truth_jets;

      // Finally, set the hardest truth jet.

      set_hardest_truth_jet();
   }
   else if (data_source() == 2) {   // 2 => MC_RECO.

      // Recluster all reco particles to get reco jets.
      
      vector<PseudoJet> reco_particles_pseudojets = convert_to_pseudojets(_mc_reco_particles);
      ClusterSequence cs(reco_particles_pseudojets, jet_def);
      vector<PseudoJet> reco_ak5_jets = sorted_by_pt(cs.inclusive_jets(0.0));

      // Create a vector of MOD::MCCalibratedJet.
      vector<MOD::MCCalibratedJet> reco_jets;
      for (unsigned i = 0; i < reco_ak5_jets.size(); i++) {
         reco_jets.push_back(MOD::MCCalibratedJet( reco_ak5_jets[i], "ak5" ));
      }

      _mc_reco_jets = reco_jets;      

      // Finally, set the hardest reco jets.

      set_hardest_reco_jet();
   }


}







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
