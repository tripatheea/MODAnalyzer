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

#include <cstdio>

#include "fastjet/ClusterSequence.hh"
#include "../interface/Event.h"

using namespace std;
using namespace fastjet;

void convert_to_one_jet(MOD::Event & event_being_read, ofstream & output_file);

int main(int argc, char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();

   int number_of_events_to_process;

   if (argc <= 3) {
        std::cerr << "ERROR: You need to supply four arguments- first, path to the input data; second, path to the output file; ; third, path to save the log file to, fourth, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./filter (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }
   else if (argc == 4) {
      // Fourth argument is missing, process everything.
      number_of_events_to_process = std::numeric_limits<int>::max();
   }
   else {
      // Fourth argument gives the number of events to process.
      number_of_events_to_process = stoi(argv[4]);
   }

   ofstream log_stream(argv[3], ios::out | ios::app );
   ofstream num_stream( std::string(argv[3]) + ".num", ios::out | ios::app );

   // connect stream buffers
   std::streambuf *cerrbuf = std::cerr.rdbuf();
   std::cerr.rdbuf(log_stream.rdbuf () );
   
   ifstream data_file(argv[1]);
   ofstream output_file(argv[2], ios::out | ios::app );
   
   cout << endl << endl << "Starting convert_to_one_jet with the following given arguments: " << endl;
   cout << "Input file: " << argv[1] << endl;
   cout << "Output file: " << argv[2] << endl;
   cout << "Number of events: ";
   if(argc == 4)
      cout << "ALL" << endl << endl;
   else
      cout << number_of_events_to_process << endl << endl;

   MOD::Event event_being_read;

   int event_serial_number = 1;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {

      if( (event_serial_number % 1000) == 0 )
         cout << "Converting event number " << event_serial_number << " to 1 jet form." << endl;

      if (event_being_read.assigned_trigger_fired()) {
        convert_to_one_jet(event_being_read, output_file);
      }

      event_being_read = MOD::Event();
      event_serial_number++;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();

   cout << "Finished converting " << (event_serial_number - 1) << " events to 1 jet form in " << elapsed_seconds << " seconds!" << endl << endl;

   // Write number of events processed to the num stream so that we can sum it up later.
   num_stream << (event_serial_number - 1) << endl;

   // Restore.
   std::cout.flush ();
   std::cout.rdbuf(cerrbuf);

}


void convert_to_one_jet(MOD::Event & event_being_read, ofstream & output_file) {

   PseudoJet trigger_jet = event_being_read.trigger_jet();

   event_being_read.convert_to_one_jet();


   if (event_being_read.trigger_jet_is_matched() && (trigger_jet.user_info<MOD::InfoCalibratedJet>().jet_quality() >= 1)) {   // Jet quality level: FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3      
      output_file << event_being_read;
   }
   
}
