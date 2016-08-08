#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>
#include <iomanip>
#include <limits>
#include <chrono>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"


#include "../interface/Event.h"
#include "../interface/Property.h"

using namespace std;
using namespace fastjet;
using namespace contrib;

void analyze_event(string mc_filename, MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number);

int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();

   int number_of_events_to_process;

   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply five arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }
   else if (argc == 3) {
      // Fifth argument is missing, process everything.
      number_of_events_to_process = std::numeric_limits<int>::max();
   }
   else {
      // Fifth argument gives the number of events to process.
      number_of_events_to_process = stoi(argv[3]);
   }

   ifstream data_file(argv[1]);
   ofstream output_file(argv[2], ios::out | ios::app);


   cout << endl << endl << "Starting analysis with the following given arguments: " << endl;
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
      
      if( (event_serial_number % 5000) == 0 )
         cout << "Processing event number " << event_serial_number << endl;

      analyze_event(argv[1], event_being_read, output_file, event_serial_number);
      
      event_being_read = MOD::Event();
      event_serial_number++;

   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_event(string mc_filename, MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number) {

   vector<MOD::Property> properties;
   
   properties.push_back( MOD::Property("mc_filename", mc_filename) );
   properties.push_back( MOD::Property("serial_number", event_serial_number) );
   properties.push_back( MOD::Property("weight", event_being_read.weight()) );

 
   

   // Now that we've calculated all observables, write them out.

   string name;
   
   int padding = 35;

   if (event_serial_number == 1) {
      for (unsigned p = 0; p < properties.size(); p++) {
         if (p > 0)
            output_file << setw(padding);
         
         output_file << "; " << properties[p].name();
      }

      output_file << endl;
   }

   for (unsigned q = 0; q < properties.size(); q++) {
      if (q > 0)
         output_file << setw(padding);
      output_file <<  "; " << properties[q];
   }

   output_file << endl;
   
   
}

