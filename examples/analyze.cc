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


#include "../interface/event.h"
#include "../interface/fractional_jet_multiplicity.h"
#include "../interface/property.h"

using namespace std;
using namespace fastjet;
using namespace contrib;

void analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number,  vector<double> cone_radii, vector<double> pt_cuts);

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

   vector<double> cone_radii = {0.3, 0.5, 0.7};
   vector<double> pt_cuts = {50.0, 80.0, 110.0};

   MOD::Event event_being_read;

   int event_serial_number = 1;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {
      
      if( (event_serial_number % 100) == 0 )
         cout << "Processing event number " << event_serial_number << endl;

      analyze_event(event_being_read, output_file, event_serial_number, cone_radii, pt_cuts);
      
      event_being_read = MOD::Event();
      event_serial_number++;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, vector<double> cone_radii, vector<double> pt_cuts) {

   
   
   vector<MOD::Property> properties;

   properties.push_back(MOD::Property("# Entry", "  Entry"));

   properties.push_back(MOD::Property("Event_Number", event_being_read.event_number()));
   properties.push_back(MOD::Property("Run_Number", event_being_read.run_number()));
   
   

   MOD::CalibratedJet hardest_jet = event_being_read.hardest_jet();
   PseudoJet corrected_hardest_pseudojet = hardest_jet.pseudojet() * hardest_jet.JEC();

   properties.push_back(MOD::Property("Hardest_pT", hardest_jet.pseudojet().pt()));
   properties.push_back(MOD::Property("Corr_Hardest_pT", corrected_hardest_pseudojet.pt()));

   properties.push_back(MOD::Property("Prescale", event_being_read.assigned_trigger_prescale()));
   properties.push_back(MOD::Property("Trigger_Name", event_being_read.assigned_trigger_name()));


   // Run AK5 clustering with FastJet to get zg value.

   JetDefinition jet_def(antikt_algorithm, 0.5);
   ClusterSequence cs(event_being_read.pseudojets(), jet_def);
   vector<PseudoJet> ak5_jets = sorted_by_pt(cs.inclusive_jets());

   if (ak5_jets.size() > 0) {
      PseudoJet hardest_jet = ak5_jets[0];
      // hardest_jet *= event_being_read.hardest_jet_JEC();

      double beta = 0;

      SoftDrop soft_drop(beta, 0.05);
      PseudoJet soft_drop_jet = soft_drop(hardest_jet);
      double zg_05 = soft_drop_jet.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_05", zg_05));

      SoftDrop soft_drop_1(beta, 0.1);
      PseudoJet soft_drop_jet_1 = soft_drop_1(hardest_jet);
      double zg_1 = soft_drop_jet_1.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_1", zg_1));  

      SoftDrop soft_drop_2(beta, 0.2);
      PseudoJet soft_drop_jet_2 = soft_drop_2(hardest_jet);
      double zg_2 = soft_drop_jet_2.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_2", zg_2));  
   }
   else {
      properties.push_back(MOD::Property("zg_05", 0.0));      
      properties.push_back(MOD::Property("zg_1", 0.0));      
      properties.push_back(MOD::Property("zg_2", 0.0));      
   }


   string name;
   
   int padding = 20;

   if (event_serial_number == 1) {
      for (unsigned p = 0; p < properties.size(); p++) {
         
         if (p > 0)
            output_file << setw(padding);
         
         output_file << properties[p].name();
      }

      output_file << endl;
   }

   for (unsigned q = 0; q < properties.size(); q++) {
      if (q > 0)
         output_file << setw(padding);
      output_file << properties[q];
   }

   output_file << endl;

   
}