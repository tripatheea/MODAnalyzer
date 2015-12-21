#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>
#include <iomanip>

#include "fastjet/ClusterSequence.hh"
#include "../interface/event.h"
#include "../interface/fractional_jet_multiplicity.h"

using namespace std;
using namespace fastjet;


void write(MOD::Event & event_being_read, ofstream & output_file);


int main(int argc, char * argv[]) {
   
   int number_of_events_to_process;

   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply three arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./duplication (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
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



   cout << endl << endl << "Starting validation with the following given arguments: " << endl;
   cout << "Input file: " << argv[1] << endl;
   cout << "Output file: " << argv[2] << endl;
   cout << "Number of events: ";
   if(argc == 3)
      cout << "ALL" << endl << endl;
   else
      cout << number_of_events_to_process << endl << endl;
   
 
   MOD::Event event_being_read;

   int event_serial_number = 1;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {

      if( (event_serial_number % 100) == 0 )
         cout << "Writing out event number " << event_serial_number << endl;
     
      write(event_being_read, output_file);

      
      event_being_read = MOD::Event();
      event_serial_number++;
   }
   
   cout << "Writing jets complete!" << endl << endl;
   

}



void write(MOD::Event & event_being_read, ofstream & output_file) {
   output_file << event_being_read;
}