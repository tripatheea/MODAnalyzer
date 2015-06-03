#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>
#include <iomanip>
#include <chrono>

#include "fastjet/ClusterSequence.hh"

#include "event.cc"
#include "fractional_jet_multiplicity.cc"

using namespace std;


bool read_event(ifstream & data_file, MODEvent & event);
bool analyze_event(MODEvent & event_being_read, ofstream & output_file, vector<double> cone_radii, vector<double> pt_cuts);
vector<string> split(string const &input);

int main() {
	ifstream data_file("../data/pfcandidates.dat");

	ofstream output_file("../data/antikt_multiplicities.dat", ios::out);
	
	vector<double> cone_radii = {0.3, 0.5, 0.7};
	vector<double> pt_cuts = {50.0, 80.0, 110.0};

	MODEvent * event_being_read = new MODEvent();

	output_file << "# Event_Number     Run_Number     N_tilde     Jet_Size          Trigger_Name          Fired?     Prescale_1     Prescale_2     Cone_Radius     pT_Cut     Hardest_pT" << endl;
	
	int event_serial_number = 1;
	while(read_event(data_file, * event_being_read)) {

		// cout << event_being_read->make_string();
		analyze_event( * event_being_read, output_file, cone_radii, pt_cuts);
		
		cout << "Processing event number " << event_serial_number << endl;
		
		delete event_being_read;
		MODEvent * event_being_read = new MODEvent();

		event_serial_number++;
	}
}


vector<string> split(string const &input) { 
    istringstream buffer(input);
    vector<string> ret((istream_iterator<string>(buffer)), istream_iterator<string>());
    return ret;
}

bool read_event(ifstream & data_file, MODEvent & event_being_read) {

	string line;
	while(getline(data_file, line)) {
		istringstream iss(line);

		vector<string> components = split(line);

		if (components[0] == "BeginEvent") {
			event_being_read.set_event_number(stoi(components[4]));
			event_being_read.set_run_number(stoi(components[2]));
			event_being_read.set_particles_trigger_type("PFC");
		}
		else if (components[0] == "PFC") {
			try {
				event_being_read.add_particle(line);
			}
			catch (exception& e) {
				throw runtime_error("Invalid file format!");
				cout << "Something went wrong!" << endl;
			}
		}
		else if (components[0] == "trig") {
			try {
				event_being_read.add_trigger(line);
			}
			catch (exception& e) {
				throw runtime_error("Invalid file format!");
				cout << "Something went wrong!" << endl;
			}
		}
		else if (components[0] == "EndEvent") {
			return true;
		}
	}

	return false;
}

bool analyze_event(MODEvent & event_being_read, ofstream & output_file, vector<double> cone_radii, vector<double> pt_cuts) {

	// Retrieve the assigned trigger and store information about that trigger (prescales, fired or not).

	// Also calculate everything and record those along with the trigger information.

	string assigned_trigger_name = event_being_read.assigned_trigger_name();
	MODTrigger assigned_trigger = event_being_read.trigger_by_name(assigned_trigger_name);

	pair<int, int> prescales = assigned_trigger.prescale_pair();
	
	bool fired = assigned_trigger.fired();
	

	int prescale_1 = prescales.first;
	int prescale_2 = prescales.second;

	double hardest_pt = event_being_read.hardest_pt();

	// Calculate everything for each value of R and pt_cut.

	for(unsigned int r = 0; r < cone_radii.size(); r++) {
		for (unsigned int p = 0; p < pt_cuts.size(); p++) {

			// Calculate N_tilde.
			FractionalJetMultiplicity n_tilde_1 = FractionalJetMultiplicity(cone_radii[r], pt_cuts[p]);
			double N_tilde = n_tilde_1.calculate_n_tilde(event_being_read.pseudojets());
			
			// Calculate jet size (fastjet)
			JetDefinition jet_def(antikt_algorithm, cone_radii[r]);
			ClusterSequence cs(event_being_read.pseudojets(), jet_def);
			vector<PseudoJet> jets = cs.inclusive_jets(pt_cuts[p]);

			output_file << setw(12) << event_being_read.event_number()
						<< setw(15) << event_being_read.run_number()
						<< setw(14) << showpoint << setprecision(6) << N_tilde
						<< setw(9) << jets.size()
						<< setw(35) << assigned_trigger_name
						<< setw(5) << fired
						<< setw(12) << prescale_1
						<< setw(15) << prescale_2
						<< setw(18) << setprecision(2) << cone_radii[r]
						<< setw(12) << noshowpoint << setprecision(3) << pt_cuts[p]
						<< setw(16) << showpoint << setprecision(8) << hardest_pt
						<< endl;					
		}
	}
}