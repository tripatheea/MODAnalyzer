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
#include <cmath>

#include "fastjet/ClusterSequence.hh"

#include "fastjet/contrib/SoftDrop.hh"

#include "../interface/event.h"

using namespace std;
using namespace fastjet;
using namespace contrib;

void analyze_zg(MOD::Event & event_being_read, ofstream & output_file, vector<double> z_cuts);

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
   ofstream output_file(argv[2], ios::out | ios::app);
   
   cout << endl << endl << "Starting analysis with the following given arguments: " << endl;
   cout << "Input file: " << argv[1] << endl;
   cout << "Output file: " << argv[2] << endl;
   cout << "Number of events: ";
   if(argc == 3)
      cout << "ALL" << endl << endl;
   else
      cout << number_of_events_to_process << endl << endl;

   vector<double> z_cuts = {0.05, 0.1, 0.2};

   MOD::Event event_being_read;

   output_file << "# Entries         z_cut    hardest_jet_pt                 zg    prescale " << endl;

   int event_serial_number = 1;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {
      
      if( (event_serial_number % 100) == 0 )
         cout << "Processing event number " << event_serial_number << endl;

      analyze_zg(event_being_read, output_file, z_cuts);
      event_being_read = MOD::Event();
      event_serial_number++;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_zg(MOD::Event & event_being_read, ofstream & output_file, vector<double> z_cuts) {

   // Retrieve the assigned trigger and store information about that trigger (prescales, fired or not).
   // Also calculate everything and record those along with the trigger information.

   bool fired = event_being_read.assigned_trigger_fired();
   int prescale = event_being_read.assigned_trigger_prescale();

   // Calculate everything for each value of z_cut and pt_cut.

   if (fired) {
      for(unsigned int z = 0; z < z_cuts.size(); z++) {
            
         // Run AK5 clustering with FastJet.
         JetDefinition jet_def(antikt_algorithm, 0.5);
         ClusterSequence cs(event_being_read.pseudojets(), jet_def);
         vector<PseudoJet> ak5_jets = sorted_by_pt(cs.inclusive_jets());

         if (ak5_jets.size() > 0) {
            PseudoJet hardest_jet = ak5_jets[0];
            
            // hardest_jet *= event_being_read.hardest_jet_JEC();

            double beta = 0;

            SoftDrop soft_drop(beta, z_cuts[z]);
            PseudoJet soft_drop_jet = soft_drop(hardest_jet);
            double zg = soft_drop_jet.structure_of<SoftDrop>().symmetry();

            output_file << "   ENTRY"
                     << setw(15) << z_cuts[z]
                     << setw(18) << hardest_jet.pt()
                     << setw(19) << showpoint << setprecision(8) << zg
                     << setw(12) << prescale
                     << endl;      
         }       
      }
   }
}
