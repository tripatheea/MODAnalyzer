#include "mc_event.h"


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

const vector<MOD::PFCandidate> & MOD::Event::particles() const {
   return _particles;
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

string MOD::Event::make_string() const {
   stringstream file_to_write;
   
   file_to_write << "BeginEvent Version " << _version << " " << _data_type.first << " " << _data_type.second << endl;
   

   // First, write out conditions.
   file_to_write << _condition.make_header_string();
   file_to_write << _condition;

   

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



void MOD::Event::establish_properties() {
   
   // Cluster PFCandidates into AK5 Jets using FastJet.
   vector<MOD::PFCandidate> pfcandidates = particles();
   vector<PseudoJet> pfcandidates_pseudojets = MOD::convert_to_pseudojets(pfcandidates);

   JetDefinition jet_def(antikt_algorithm, 0.5);
   ClusterSequence cs(pfcandidates_pseudojets, jet_def);
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






vector<MOD::CalibratedJet> MOD::Event::apply_eta_cut(vector<MOD::CalibratedJet> jets, double eta_cut) const {

   vector<MOD::CalibratedJet> jets_with_eta_cut;

   for (unsigned i = 0; i < jets.size(); i++) {
      if ( abs(jets[i].uncorrected_pseudojet().eta()) < eta_cut ) {
         jets_with_eta_cut.push_back(jets[i]);
      }
   }

   return jets_with_eta_cut;
}



namespace MOD {
   ostream& operator<< (ostream& os, const Event& event) {
      os << event.make_string();
      return os;
   }
}
