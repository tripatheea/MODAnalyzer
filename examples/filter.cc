#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>
#include <iomanip>

#include "../src/event.cc"
#include "../src/fractional_jet_multiplicity.cc"

using namespace std;

bool filter_events(MODEvent & event_being_read, ofstream & output_file);

int main(int argc, char * argv[]) {
	
	if (argc != 2) {
        std::cerr << "ERROR: You need to supply path to the input data. The path has to be either absolute or relative to the bin directory." << std::endl;
        return 1;
    }

	ifstream data_file(argv[1]);
	ofstream output_file("../data/filtered_events.dat", ios::out);
	
	MODEvent * event_being_read = new MODEvent();

	int event_serial_number = 1;
	while(event_being_read->read_event(data_file)) {

		filter_events( * event_being_read, output_file);
		
		cout << "Filtering event number " << event_serial_number << endl;
		
		delete event_being_read;
		MODEvent * event_being_read = new MODEvent();

		event_serial_number++;
	}
}


bool filter_events(MODEvent & event_being_read, ofstream & output_file) {

	// Retrieve the assigned trigger and store information about that trigger (prescales, fired or not).

	// Also calculate everything and record those along with the trigger information.

	string assigned_trigger_name = event_being_read.assigned_trigger_name();
	const MODTrigger assigned_trigger = event_being_read.trigger_by_name(assigned_trigger_name);
	
	if(assigned_trigger.fired()) {
		output_file << event_being_read.make_string();
	}
	
}