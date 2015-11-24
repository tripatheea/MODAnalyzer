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

void analyze_qcd_beta(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, std::string);

double calculate_rho(double R, double m, double pT);


int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();

   int number_of_events_to_process;

   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply three arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }
   else if (argc == 4) {
      // Third argument is missing, process everything.
      number_of_events_to_process = std::numeric_limits<int>::max();
   }
   else {
      // Third argument gives the number of events to process.
      number_of_events_to_process = stoi(argv[4]);
   }

   ifstream data_file(argv[1]);
   ofstream output_file(argv[2], ios::out | ios::app);
   std::string reco_or_truth = argv[3];

   
   cout << endl << endl << "Starting analysis with the following given arguments: " << endl;
   cout << "Input file: " << argv[1] << endl;
   cout << "Output file: " << argv[2] << endl;
   cout << "Reco vs. Truth: " << argv[3] << endl;
   cout << "Number of events: ";

   if(argc == 4)
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

      
      // Write out version info in the output file for the "syncing plots" thing to work (as it needs to figure out which directory to put things into).
      if (event_serial_number == 1)
         output_file << "%" << " Version " << event_being_read.version() << endl;

      analyze_qcd_beta(event_being_read, output_file, event_serial_number, reco_or_truth);
      
      
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



void analyze_qcd_beta(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, std::string reco_or_truth) {

   JetDefinition jet_def_cambridge(cambridge_algorithm, jet_def_cambridge.max_allowable_R);

   vector<fastjet::PseudoJet> hardest_jet_constituents;
   if (reco_or_truth == "reco") {
      hardest_jet_constituents = event_being_read.hardest_mc_reco_jet_constituents();   
   }
   else {
      hardest_jet_constituents = event_being_read.hardest_mc_truth_jet_constituents();     
   }
   
   

   ClusterSequence cs(hardest_jet_constituents, jet_def_cambridge);
   if (cs.inclusive_jets().size() == 0) {
      return;
   }

   fastjet::PseudoJet hardest_jet = cs.inclusive_jets()[0];

   vector<MOD::Property> properties;

   properties.push_back(MOD::Property("# Entry", "  Entry"));
   properties.push_back(MOD::Property("Hardest_pT", hardest_jet.pt()));
   
   
   
   double beta = 0;
   
   SoftDrop soft_drop_10(beta, 0.10);
   PseudoJet soft_drop_jet_10 = soft_drop_10(hardest_jet);
   properties.push_back(MOD::Property("Rg_10", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_10", soft_drop_jet_10.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_11(beta, 0.11);
   PseudoJet soft_drop_jet_11 = soft_drop_11(hardest_jet);
   properties.push_back(MOD::Property("Rg_11", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_11", soft_drop_jet_11.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_12(beta, 0.12);
   PseudoJet soft_drop_jet_12 = soft_drop_12(hardest_jet);
   properties.push_back(MOD::Property("Rg_12", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_12", soft_drop_jet_12.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_13(beta, 0.13);
   PseudoJet soft_drop_jet_13 = soft_drop_13(hardest_jet);
   properties.push_back(MOD::Property("Rg_13", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_13", soft_drop_jet_13.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_14(beta, 0.14);
   PseudoJet soft_drop_jet_14 = soft_drop_14(hardest_jet);
   properties.push_back(MOD::Property("Rg_14", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_14", soft_drop_jet_14.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_15(beta, 0.15);
   PseudoJet soft_drop_jet_15 = soft_drop_15(hardest_jet);
   properties.push_back(MOD::Property("Rg_15", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_15", soft_drop_jet_15.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_16(beta, 0.16);
   PseudoJet soft_drop_jet_16 = soft_drop_16(hardest_jet);
   properties.push_back(MOD::Property("Rg_16", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_16", soft_drop_jet_16.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_17(beta, 0.17);
   PseudoJet soft_drop_jet_17 = soft_drop_17(hardest_jet);
   properties.push_back(MOD::Property("Rg_17", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_17", soft_drop_jet_17.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_18(beta, 0.18);
   PseudoJet soft_drop_jet_18 = soft_drop_18(hardest_jet);
   properties.push_back(MOD::Property("Rg_18", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_18", soft_drop_jet_18.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_19(beta, 0.19);
   PseudoJet soft_drop_jet_19 = soft_drop_19(hardest_jet);
   properties.push_back(MOD::Property("Rg_19", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_19", soft_drop_jet_19.structure_of<SoftDrop>().symmetry()));

   SoftDrop soft_drop_20(beta, 0.20);
   PseudoJet soft_drop_jet_20 = soft_drop_20(hardest_jet);
   properties.push_back(MOD::Property("Rg_20", soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
   properties.push_back(MOD::Property("zg_20", soft_drop_jet_20.structure_of<SoftDrop>().symmetry()));


   // Just get "charged" particles.

   vector<fastjet::PseudoJet> charged_constituents = MOD::filter_charged(hardest_jet_constituents);
   ClusterSequence cs_charged(charged_constituents, jet_def_cambridge);

   if (cs_charged.inclusive_jets().size() == 0) {

      properties.push_back(MOD::Property("chrg_Rg_10", -1.));
      properties.push_back(MOD::Property("chrg_zg_10", -1.));

      properties.push_back(MOD::Property("chrg_Rg_11", -1.));
      properties.push_back(MOD::Property("chrg_zg_11", -1.));
      
      properties.push_back(MOD::Property("chrg_Rg_12", -1.));
      properties.push_back(MOD::Property("chrg_zg_12", -1.));
      
      properties.push_back(MOD::Property("chrg_Rg_13", -1.));
      properties.push_back(MOD::Property("chrg_zg_13", -1.));
      
      properties.push_back(MOD::Property("chrg_Rg_14", -1.));
      properties.push_back(MOD::Property("chrg_zg_14", -1.));
      
      properties.push_back(MOD::Property("chrg_Rg_15", -1.));
      properties.push_back(MOD::Property("chrg_zg_15", -1.));
      
      properties.push_back(MOD::Property("chrg_Rg_16", -1.));
      properties.push_back(MOD::Property("chrg_zg_16", -1.));
      
      properties.push_back(MOD::Property("chrg_Rg_17", -1.));
      properties.push_back(MOD::Property("chrg_zg_17", -1.));
      
      properties.push_back(MOD::Property("chrg_Rg_18", -1.));
      properties.push_back(MOD::Property("chrg_zg_18", -1.));
      
      properties.push_back(MOD::Property("chrg_Rg_19", -1.));
      properties.push_back(MOD::Property("chrg_zg_19", -1.));
      
      properties.push_back(MOD::Property("chrg_Rg_20", -1.));
      properties.push_back(MOD::Property("chrg_zg_20", -1.));
   }
   else {

      PseudoJet hardest_charged_jet = cs.inclusive_jets()[0];
      
      SoftDrop charged_soft_drop_10(beta, 0.10);
      PseudoJet charged_soft_drop_jet_10 = charged_soft_drop_10(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_10", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_10", charged_soft_drop_jet_10.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_11(beta, 0.11);
      PseudoJet charged_soft_drop_jet_11 = charged_soft_drop_11(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_11", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_11", charged_soft_drop_jet_11.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_12(beta, 0.12);
      PseudoJet charged_soft_drop_jet_12 = charged_soft_drop_12(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_12", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_12", charged_soft_drop_jet_12.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_13(beta, 0.13);
      PseudoJet charged_soft_drop_jet_13 = charged_soft_drop_13(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_13", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_13", charged_soft_drop_jet_13.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_14(beta, 0.14);
      PseudoJet charged_soft_drop_jet_14 = charged_soft_drop_14(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_14", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_14", charged_soft_drop_jet_14.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_15(beta, 0.15);
      PseudoJet charged_soft_drop_jet_15 = charged_soft_drop_15(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_15", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_15", charged_soft_drop_jet_15.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_16(beta, 0.16);
      PseudoJet charged_soft_drop_jet_16 = charged_soft_drop_16(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_16", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_16", charged_soft_drop_jet_16.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_17(beta, 0.17);
      PseudoJet charged_soft_drop_jet_17 = charged_soft_drop_17(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_17", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_17", charged_soft_drop_jet_17.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_18(beta, 0.18);
      PseudoJet charged_soft_drop_jet_18 = charged_soft_drop_18(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_18", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_18", charged_soft_drop_jet_18.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_19(beta, 0.19);
      PseudoJet charged_soft_drop_jet_19 = charged_soft_drop_19(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_19", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_19", charged_soft_drop_jet_19.structure_of<SoftDrop>().symmetry()));

      SoftDrop charged_soft_drop_20(beta, 0.20);
      PseudoJet charged_soft_drop_jet_20 = charged_soft_drop_20(hardest_charged_jet);
      properties.push_back(MOD::Property("chrg_Rg_20", charged_soft_drop_jet_10.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("chrg_zg_20", charged_soft_drop_jet_20.structure_of<SoftDrop>().symmetry()));


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