#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>

#include "fastjet/ClusterSequence.hh"
#include "event.cc"




using namespace std;

vector<string> split(string const &input);

int main() {
	ifstream data_file("pfcandidates");

	ofstream multiplicities_output_file("antikt_multiplicities.dat", ios::out);
	
	double px, py, pz, energy;	
	double R = 0.5;
	double pt_cut = 50.00;

	int run_number, event_number;

	Event * event_being_read = new Event();
	
	int event_serial_number = 1;

	string line;
	while(getline(data_file, line)) {
		istringstream iss(line);

		cout << "Processing event number: " << event_serial_number << endl;

		vector<string> components = split(line);

		if (components[0] == "BeginEvent") {
			run_number = stoi(components[2]);
			event_number = stoi(components[4]);

			event_being_read->assign_event_number(event_number);
			event_being_read->assign_run_number(run_number);
			event_being_read->assign_particles_type("PFC");
		}
		else if (components[0] == "PFC") {
			try {
				// cout << components[1] << endl;
				cout.precision(8);
				event_being_read->add_particle(stod(components[1]), stod(components[2]), stod(components[3]), stod(components[4]), stod(components[5]), stod(components[6]));
			}
			catch (exception& e) {
				throw runtime_error("Invalid file format!");
				cout << "Something went wrong!" << endl;
			}
		}
		else if (components[0] == "EndEvent") {
			// Current event has ended. We can process it now.
			// event_being_read->write_to_file("Test.dat");

			// Now retrieve the assigned trigger and use that to determine:
			// a) whether to record the current event or not.
			// b) what prescale to use.

			Trigger assigned_trigger = event_being_read->assigned_trigger();

			if (assigned_trigger.is_valid()) {

				// Record things only if at least one trigger was fired.
				
				pair<int, int> prescales = assigned_trigger.prescales();

				int prescale_1 = prescales.first;
				int prescale_2 = prescales.second;

				// Calculate N_tilde.
				double N_tilde = event_being_read->calculate_N_tilde(R, pt_cut);
				// cout << N_tilde << endl;
				
				// Calculate jet size (fastjet)
				JetDefinition jet_def(antikt_algorithm, R);
				vector<PseudoJet> jets = event_being_read->jets(jet_def, pt_cut);


				multiplicities_output_file << N_tilde << " " << jets.size() << " " << prescale_1 << " " << prescale_2 << " " << assigned_trigger.name() << endl;
				// cout << N_tilde << " " << jets.size() << " " << prescale_1 << " " << prescale_2 << " " << assigned_trigger.name() << endl;
			}

			// We've calculated all we need. Now remove the old pointer (for the previous instance of event) and create a new one.

			delete event_being_read;
			Event * event_being_read = new Event();

			// Increment the event serial number- this is just for outputting the current event number to the user. This number isn't used internally.
			event_serial_number++;
		}
		else if (components[0] == "trig") {
			try {
				event_being_read->add_trigger(components[1], stoi(components[2]), stoi(components[3]), stoi(components[4]) == 1);
			}
			catch (exception& e) {
				throw runtime_error("Invalid file format!");
				cout << "Something went wrong!" << endl;
			}
		}

	}

	cout << "I'm done here. I'm just going to leave now." << endl;


	return 0;
}


vector<string> split(string const &input) { 
    istringstream buffer(input);
    vector<string> ret((istream_iterator<string>(buffer)), istream_iterator<string>());
    return ret;
}















// int main() {
// 	// ifstream pfCandidatesFile("../pfcandidates");
// 	ifstream pfCandidatesFile("../minBias");

// 	ofstream fmatch("antikt_multiplicities.csv", ios::out);
// 	CSVRow row;

// 	double px, py, pz, energy;	
// 	double R = 0.5;
// 	double pt_cut = 50.00;

// 	int run_number, event_number_of_next_event;
// 	int event_count = 1;
	
// 	int first_event_number = 138867139;
// 	int first_run_number = 146511;

// 	int event_number_of_event_being_processed = first_event_number;
// 	Event * current_event = new Event(first_run_number, first_event_number);
// 	int no_of_triggers = 1; 	// Depends on the analysis we're doing.



// 	while(pfCandidatesFile >> row) {

// 		run_number = stoi(row[0]);
// 		event_number_of_next_event = stoi(row[1]);
// 		px = stod(row[2]);
// 		py = stod(row[3]);
// 		pz = stod(row[4]);
// 		energy = stod(row[5]);

// 		PseudoJet current_particle = PseudoJet(px, py, pz, energy);

// 		if (event_number_of_event_being_processed != event_number_of_next_event) {

// 			cout << "Processing event number: " << event_number_of_event_being_processed << " which is # " << event_count << endl;

			
// 			// We've moved on to a different event.
// 			// That means, the class current_event contains all the particles that it's supposed to.
			
// 			// Add all the trigger info for this event.
			
// 			int line_start = (event_count - 1) * no_of_triggers + 1;

// 			vector<Trigger> current_triggers = get_trigger_info(event_number_of_event_being_processed, line_start, line_start + no_of_triggers - 1);
// 			current_event->add_triggers(current_triggers);




// 			// All triggers stored.

// 			// Now retrieve the assigned trigger and use that to determine:
// 			// a) whether to record the current event or not.
// 			// b) what prescale to use.		

// 			Trigger assigned_trigger = current_event->get_assigned_trigger();
// 			// cout << assigned_trigger.get_name() << endl;

// 			if (assigned_trigger.is_valid()) {

// 				// Record things only if at least one trigger was fired.
				
// 				pair<int, int> prescales = assigned_trigger.get_prescales();

// 				int prescale_1 = get<0>(prescales);
// 				int prescale_2 = get<1>(prescales);

// 				// Calculate N_tilde.
// 				double N_tilde = current_event->calculate_N_tilde(event_number_of_event_being_processed, R, pt_cut);
// 				// cout << N_tilde << endl;
				
// 				// Calculate jet size (fastjet)
// 				JetDefinition jet_def(antikt_algorithm, R);
// 				vector<PseudoJet> jets = current_event->get_jets(event_number_of_event_being_processed, jet_def, pt_cut);


// 				// fmatch << N_tilde << " " << jets.size() << " " << prescale_1 << " " << prescale_2 << " " << assigned_trigger.get_name() << endl;
// 				fmatch << current_event->get_hardest_pt() << "," << prescale_1 << "," << prescale_2 << " " << assigned_trigger.get_name() << endl;
// 			}

// 			// We've calculated all we need. Now delete the old pointer (for the previous instance of event) and create a new one.
// 			delete current_event;

// 			Event * current_event = new Event(run_number, event_number_of_next_event);

// 			event_number_of_event_being_processed = event_number_of_next_event;

// 			event_count++;
// 		}

// 		current_event->add_particle(px, py, pz, energy);
// 	}

// 	// Need to do all these stuff one last time for the final event.

// 	cout << "Processing event number: " << event_number_of_event_being_processed << " which is # " << event_count << ", the last event." << endl;

// 	// Add all the trigger info for this event.
			
// 	int line_start = (event_count - 1) * no_of_triggers + 1;

// 	current_event->add_triggers(get_trigger_info(event_number_of_event_being_processed, line_start, line_start + no_of_triggers - 1));

// 	// All triggers stored.

// 	// Now retrieve the assigned trigger and use that to determine:
// 	// a) whether to record the current event or not.
// 	// b) what prescale to use.

// 	Trigger assigned_trigger = current_event->get_assigned_trigger();

// 	if (assigned_trigger.is_valid()) {

// 		// Record things only if at least one trigger was fired.

// 		string name = assigned_trigger.get_name();

// 		pair<int, int> prescales = assigned_trigger.get_prescales();

// 		int prescale_1 = get<0>(prescales);
// 		int prescale_2 = get<1>(prescales);


// 		// Calculate N_tilde.
// 		double N_tilde = current_event->calculate_N_tilde(event_number_of_event_being_processed, R, pt_cut);

// 		// Calculate jet size (fastjet)
// 		JetDefinition jet_def(antikt_algorithm, R);
// 		vector<PseudoJet> jets = current_event->get_jets(event_number_of_event_being_processed, jet_def, pt_cut);
		
// 		// fmatch <<  N_tilde << "," << jets.size() << "," << prescale_1 << "," << prescale_2 << " " << assigned_trigger.get_name() << endl;
// 		fmatch << current_event->get_hardest_pt() << "," << prescale_1 << "," << prescale_2 << " " << assigned_trigger.get_name() << endl;
// 	}

// }