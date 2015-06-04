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

bool filter_events(MOD::Event & event_being_read, ofstream & output_file);

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
   
   MOD::Event event_being_read;

   int event_serial_number = 1;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {

      cout << "Filtering event number " << event_serial_number << endl;

      filter_events(event_being_read, output_file);
      event_being_read = Event();
      event_serial_number++;
   }
}


bool filter_events(MOD::Event & event_being_read, ofstream & output_file) {
   if(event_being_read.assigned_trigger_fired()) {
      output_file << event_being_read.make_string();
   }
}