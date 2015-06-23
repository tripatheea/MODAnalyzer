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

#include "../interface/event.h"
#include "../interface/fractional_jet_multiplicity.h"

using namespace std;
using namespace fastjet;

void corrected_ak5_spectrum(MOD::Event & event_being_read, ofstream & output_file);
<<<<<<< HEAD
bool pseudojets_compare(PseudoJet a, PseudoJet b);
=======
>>>>>>> fb4fa903be032715b1351db53399a0357c1f7cf4

int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();

   int number_of_events_to_process;

   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply three arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
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
<<<<<<< HEAD
   ofstream output_file(argv[2], ios::out | ios::app);
=======
   ofstream output_file(argv[2], ios::out);
>>>>>>> fb4fa903be032715b1351db53399a0357c1f7cf4
   
   cout << endl << endl << "Starting analysis with the following given arguments: " << endl;
   cout << "Input file: " << argv[1] << endl;
   cout << "Output file: " << argv[2] << endl;
   cout << "Number of events: ";
   if(argc == 3)
      cout << "ALL" << endl << endl;
   else
      cout << number_of_events_to_process << endl << endl;

   vector<double> cone_radii = {0.3, 0.5, 0.7};
   vector<double> pt_cuts = {50.0, 80.0, 110.0};

   MOD::Event event_being_read;

<<<<<<< HEAD
   output_file << "# Entries  uncorrected_pt  corrected_pt  prescale" << endl;
=======
   output_file << "# Entries  uncorrected_pt  corrected_pt  fired  prescale" << endl;
>>>>>>> fb4fa903be032715b1351db53399a0357c1f7cf4

   int event_serial_number = 1;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {
      
      if( (event_serial_number % 100) == 0 )
         cout << "Processing event number " << event_serial_number << endl;

      corrected_ak5_spectrum(event_being_read, output_file);
      event_being_read = MOD::Event();
      event_serial_number++;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void corrected_ak5_spectrum(MOD::Event & event_being_read, ofstream & output_file) {

   // Retrieve the assigned trigger and store information about that trigger (prescales, fired or not).
   // Also calculate everything and record those along with the trigger information.
   
   string assigned_trigger_name = event_being_read.assigned_trigger_name();
   bool fired = event_being_read.assigned_trigger_fired();
   int prescale = event_being_read.assigned_trigger_prescale();
   
<<<<<<< HEAD
   vector<fastjet::PseudoJet> corrected_jets = event_being_read.corrected_calibrated_pseudojets_ak5();
   vector<fastjet::PseudoJet> jets = event_being_read.calibrated_pseudojets_ak5();

   if ((fired) && (jets.size() > 0)) {
      // Sort both jets.
      sort(corrected_jets.begin(), corrected_jets.end(), pseudojets_compare);
      sort(jets.begin(), jets.end(), pseudojets_compare);

      
      output_file << "   ENTRY"
                  << setw(17) << jets[0].pt()
                  << setw(14) << corrected_jets[0].pt()
                  << setw(10) << prescale
                  << endl;
   }
}

bool pseudojets_compare(PseudoJet a, PseudoJet b) {
   if (a.pt() > b.pt())
      return true;
   return false;
=======

   vector<fastjet::PseudoJet> corrected_jets = event_being_read.corrected_calibrated_pseudojets_ak5();
   vector<fastjet::PseudoJet> jets = event_being_read.calibrated_pseudojets_ak5();

   for (unsigned i = 0; i < corrected_jets.size(); i++) {
      output_file << "   ENTRY"
                  << setw(17) << jets[i].pt()
                  << setw(14) << corrected_jets[i].pt()
                  << setw(7) << fired
                  << setw(10) << prescale
                  << endl;             
   }
   
>>>>>>> fb4fa903be032715b1351db53399a0357c1f7cf4
}