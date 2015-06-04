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

bool duplicate_events(MODEvent & event_being_read, ofstream & output_file);

int main(int argc, char * argv[]) {
   
   int number_of_events_to_process;

   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply three arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory." << std::endl;
        return 1;
   }
   else if (argc == 3) {
      // Third argument is missing, process everything.
      number_of_events_to_process = std::numeric_limits<int>::max();
   }
   else {
      // Third argument gives the number of events to process.
      number_of_events_to_process = stoi(argv[3]);
   }

   ifstream data_file(argv[1]);
   ofstream output_file(argv[2], ios::out);
   
   MODEvent event_being_read;

   int event_serial_number = 1;
   while(event_being_read.read_event(data_file)) {

      cout << "Duplicating event number " << event_serial_number << endl;

      duplicate_events(event_being_read, output_file);
      event_being_read = MODEvent();
      event_serial_number++;
   }
}


bool duplicate_events(MODEvent & event_being_read, ofstream & output_file) {

   // Retrieve the assigned trigger and store information about that trigger (prescales, fired or not).

   // Also calculate everything and record those along with the trigger information.

   string assigned_trigger_name = event_being_read.assigned_trigger_name();
   const MODTrigger assigned_trigger = event_being_read.trigger_by_name(assigned_trigger_name);
   
   output_file << event_being_read;
   
}