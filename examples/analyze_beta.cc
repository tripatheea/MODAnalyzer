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

void analyze_qcd_beta(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number,  vector<double> cone_radii, vector<double> pt_cuts);

double calculate_rho(double R, double m, double pT);


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

      if (event_being_read.jet_quality_cut("loose") && abs(event_being_read.hardest_corrected_jet().pseudojet().eta()) < 2.4) {
         analyze_qcd_beta(event_being_read, output_file, event_serial_number, cone_radii, pt_cuts);
      }
      else {
         // cout << "I reject this one!" << endl;
      }
      
      event_being_read = MOD::Event();
      event_serial_number++;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}






double calculate_rho(double R, double m, double pT) {
   return m * m / (pT * pT * R * R);
}



void analyze_qcd_beta(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, vector<double> cone_radii, vector<double> pt_cuts) {



   vector<MOD::Property> properties;

   properties.push_back(MOD::Property("# Entry", "  Entry"));

   MOD::CalibratedJet hardest_corrected_jet = event_being_read.hardest_corrected_jet();

   properties.push_back(MOD::Property("Cor_Hardest_pT", hardest_corrected_jet.pseudojet().pt()));

   properties.push_back(MOD::Property("Prescale", event_being_read.assigned_trigger_prescale()));


   // Run AK5 clustering with FastJet to get zg value.

   JetDefinition jet_def(antikt_algorithm, 0.5);
   ClusterSequence cs(event_being_read.pseudojets(), jet_def);
   vector<PseudoJet> ak5_jets = sorted_by_pt(cs.inclusive_jets(3.0));

   if ((unsigned) ak5_jets.size() > 0) {
      PseudoJet hardest_jet = ak5_jets[0];

      double beta = 0;

      
      SoftDrop soft_drop_10(beta, 0.10);
      PseudoJet soft_drop_jet_10 = soft_drop_10(hardest_jet);
      properties.push_back(MOD::Property("rho_10", calculate_rho(0.5, soft_drop_jet_10.m(), soft_drop_jet_10.pt())));
      properties.push_back(MOD::Property("zg_10", soft_drop_jet_10.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_11(beta, 0.11);
      PseudoJet soft_drop_jet_11 = soft_drop_11(hardest_jet);
      properties.push_back(MOD::Property("rho_11", calculate_rho(0.5, soft_drop_jet_11.m(), soft_drop_jet_11.pt())));
      properties.push_back(MOD::Property("zg_11", soft_drop_jet_11.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_12(beta, 0.12);
      PseudoJet soft_drop_jet_12 = soft_drop_12(hardest_jet);
      properties.push_back(MOD::Property("rho_12", calculate_rho(0.5, soft_drop_jet_12.m(), soft_drop_jet_12.pt())));
      properties.push_back(MOD::Property("zg_12", soft_drop_jet_12.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_13(beta, 0.13);
      PseudoJet soft_drop_jet_13 = soft_drop_13(hardest_jet);
      properties.push_back(MOD::Property("rho_13", calculate_rho(0.5, soft_drop_jet_13.m(), soft_drop_jet_13.pt())));
      properties.push_back(MOD::Property("zg_13", soft_drop_jet_13.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_14(beta, 0.14);
      PseudoJet soft_drop_jet_14 = soft_drop_14(hardest_jet);
      properties.push_back(MOD::Property("rho_14", calculate_rho(0.5, soft_drop_jet_14.m(), soft_drop_jet_14.pt())));
      properties.push_back(MOD::Property("zg_14", soft_drop_jet_14.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_15(beta, 0.15);
      PseudoJet soft_drop_jet_15 = soft_drop_15(hardest_jet);
      properties.push_back(MOD::Property("rho_15", calculate_rho(0.5, soft_drop_jet_15.m(), soft_drop_jet_15.pt())));
      properties.push_back(MOD::Property("zg_15", soft_drop_jet_15.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_16(beta, 0.16);
      PseudoJet soft_drop_jet_16 = soft_drop_16(hardest_jet);
      properties.push_back(MOD::Property("rho_16", calculate_rho(0.5, soft_drop_jet_16.m(), soft_drop_jet_16.pt())));
      properties.push_back(MOD::Property("zg_16", soft_drop_jet_16.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_17(beta, 0.17);
      PseudoJet soft_drop_jet_17 = soft_drop_17(hardest_jet);
      properties.push_back(MOD::Property("rho_17", calculate_rho(0.5, soft_drop_jet_17.m(), soft_drop_jet_17.pt())));
      properties.push_back(MOD::Property("zg_17", soft_drop_jet_17.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_18(beta, 0.18);
      PseudoJet soft_drop_jet_18 = soft_drop_18(hardest_jet);
      properties.push_back(MOD::Property("rho_18", calculate_rho(0.5, soft_drop_jet_18.m(), soft_drop_jet_18.pt())));
      properties.push_back(MOD::Property("zg_18", soft_drop_jet_18.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_19(beta, 0.19);
      PseudoJet soft_drop_jet_19 = soft_drop_19(hardest_jet);
      properties.push_back(MOD::Property("rho_19", calculate_rho(0.5, soft_drop_jet_19.m(), soft_drop_jet_19.pt())));
      properties.push_back(MOD::Property("zg_19", soft_drop_jet_19.structure_of<SoftDrop>().symmetry()));

      SoftDrop soft_drop_20(beta, 0.20);
      PseudoJet soft_drop_jet_20 = soft_drop_20(hardest_jet);
      properties.push_back(MOD::Property("rho_20", calculate_rho(0.5, soft_drop_jet_20.m(), soft_drop_jet_20.pt())));
      properties.push_back(MOD::Property("zg_20", soft_drop_jet_20.structure_of<SoftDrop>().symmetry()));

      
   }
  
   

   string name;
   
   int padding = 30;

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