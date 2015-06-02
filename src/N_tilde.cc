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

using namespace std;


bool read_event(ifstream & data_file, MODEvent & event);
bool analyze_event(MODEvent & event_being_read, ofstream & output_file, vector<double> cone_radii, vector<int> pt_cuts);
vector<string> split(string const &input);

int main() {
	ifstream data_file("pfcandidates.dat");

	ofstream output_file("antikt_multiplicities.dat", ios::out);
	
	vector<double> cone_radii = {0.3, 0.5, 0.7};
	vector<int> pt_cuts = {50, 80, 110};

	MODEvent * event_being_read = new MODEvent();

	output_file << "# Event_Number     Run_Number     N_tilde     Jet_Size          Trigger_Name          Fired?     Prescale_1     Prescale_2     Cone_Radius     pT_Cut     Hardest_pT" << endl;
	
	int event_serial_number = 1;
	while(read_event(data_file, * event_being_read)) {
		// event_being_read.write_to_file("Test.dat");
		analyze_event( * event_being_read, output_file, cone_radii, pt_cuts);
		
		// cout << "Processing event number " << event_serial_number << endl;
		
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

			int run_number = stoi(components[2]);
			int event_number = stoi(components[4]);

			event_being_read.assign_event_number(event_number);
			event_being_read.assign_run_number(run_number);
			event_being_read.assign_particles_type("PFC");
		}
		else if (components[0] == "PFC") {
			try {
				event_being_read.add_particle(stod(components[1]), stod(components[2]), stod(components[3]), stod(components[4]), stod(components[5]), stoi(components[6]));
			}
			catch (exception& e) {
				throw runtime_error("Invalid file format!");
				cout << "Something went wrong!" << endl;
			}
		}
		else if (components[0] == "trig") {
			try {
				event_being_read.add_trigger(components[1], stoi(components[2]), stoi(components[3]), stoi(components[4]) == 1);
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

bool analyze_event(MODEvent & event_being_read, ofstream & output_file, vector<double> cone_radii, vector<int> pt_cuts) {

	// Retrieve the assigned trigger and store information about that trigger (prescales, fired or not).

	// Also calculate everything and record those along with the trigger information.

	string assigned_trigger_name = event_being_read.assigned_trigger_name();
	MODTrigger assigned_trigger = event_being_read.trigger_by_name(assigned_trigger_name);


	pair<int, int> prescales = assigned_trigger.prescales();
	// bool fired = (assigned_trigger.is_valid()) ? assigned_trigger.fired() : 0;
	
	bool fired = assigned_trigger.fired();
	// if (assigned_trigger.fired()) {
	// 	fired = true;
	// 	cout << "Gotchha! " << fired << endl;
	// }
	

	int prescale_1 = prescales.first;
	int prescale_2 = prescales.second;

	double hardest_pt = event_being_read.hardest_pt();

	// Calculate everything for each value of R and pt_cut.

	for(unsigned int r = 0; r < cone_radii.size(); r++) {
		for (unsigned int p = 0; p < pt_cuts.size(); p++) {

			// Calculate N_tilde.
			double N_tilde = event_being_read.calculate_N_tilde(cone_radii[r], pt_cuts[p]);
			
			// Calculate jet size (fastjet)
			JetDefinition jet_def(antikt_algorithm, cone_radii[r]);
			vector<PseudoJet> jets = event_being_read.jets(jet_def, pt_cuts[p]);


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