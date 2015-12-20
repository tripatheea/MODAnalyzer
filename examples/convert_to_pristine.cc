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
#include "../interface/event.h"
#include "../interface/fractional_jet_multiplicity.h"

using namespace std;
using namespace fastjet;

void convert_to_pristine(MOD::Event & event_being_read, ofstream & output_file);
bool pseudojets_compare(PseudoJet a, PseudoJet b);

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
   
   cout << endl << endl << "Starting convert_to_pristineming with the following given arguments: " << endl;
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
         cout << "Converting event number " << event_serial_number << " to prisine form." << endl;

      if (event_being_read.assigned_trigger_fired())
         convert_to_pristine(event_being_read, output_file);

      event_being_read = MOD::Event();
      event_serial_number++;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();

   cout << "Finished converting " << (event_serial_number - 1) << " events to pristine form in " << elapsed_seconds << " seconds!" << endl << endl;

   // Write number of events processed to the num stream so that we can sum it up later.
   num_stream << (event_serial_number - 1) << endl;

   // Restore.
   std::cout.flush ();
   std::cout.rdbuf(cerrbuf);

}


void convert_to_pristine(MOD::Event & event_being_read, ofstream & output_file) {
   

   MOD::CalibratedJet trigger_jet = event_being_read.trigger_jet();
   PseudoJet closest_fastjet_jet_to_trigger_jet = event_being_read.closest_fastjet_jet_to_trigger_jet();


   PseudoJet jec_corrected_jet = closest_fastjet_jet_to_trigger_jet * trigger_jet.JEC();
   vector<PseudoJet> jec_corrected_jet_constituents = closest_fastjet_jet_to_trigger_jet.constituents();

   event_being_read.set_prescale(event_being_read.assigned_trigger_prescale());


   event_being_read.set_pristine_form(true);

   if (event_being_read.trigger_jet_is_matched() && (trigger_jet.jet_quality() >= 1)) {   // Jet quality level: FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3      
      output_file << event_being_read;
   }



   // output_file << "BeginEvent Version " << event_being_read.version() << " " << event_being_read.data_type().first << " " << event_being_read.data_type().second << " Prescale " << event_being_read.prescale() << endl;
   
   // output_file << "# PDPFC" << "              px              py              pz          energy" << endl;

   // for (unsigned i = 0; i < jec_corrected_jet_constituents.size(); i++) {
   //    output_file << "  PDPFC"
   //               << setw(16) << fixed << setprecision(8) << jec_corrected_jet_constituents[i].px()
   //               << setw(16) << fixed << setprecision(8) << jec_corrected_jet_constituents[i].py()
   //               << setw(16) << fixed << setprecision(8) << jec_corrected_jet_constituents[i].pz()
   //               << setw(16) << fixed << setprecision(8) << jec_corrected_jet_constituents[i].E()
   //               << endl;
   // }

   

   
}


bool pseudojets_compare(PseudoJet a, PseudoJet b) {
   if (a.pt() > b.pt())
      return true;
   return false;
}