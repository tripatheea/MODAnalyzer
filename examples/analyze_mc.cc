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

void mc_analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, std::string);

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
   ofstream output_file(argv[2], ios::out);
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
      
      if( (event_serial_number % 100) == 0 ) {
         cout << "Processing event number " << event_serial_number << endl;
      }

      // Write out version info in the output file for the "syncing plots" thing to work (as it needs to figure out which directory to put things into).
      if (event_serial_number == 1)
         output_file << "%" << " Version " << event_being_read.version() << endl;

      
      mc_analyze_event(event_being_read, output_file, event_serial_number, reco_or_truth);
      
      event_being_read = MOD::Event();
      event_serial_number++;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void mc_analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, std::string reco_or_truth) {


   JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);

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
   // properties.push_back(MOD::Property("Hardest_pT", event_being_read.hardest_mc_truth_jet().pseudojet().pt()));
   

   // Get all 


   SoftDrop soft_drop(0.0, 0.05);
   PseudoJet soft_drop_jet = soft_drop(hardest_jet);
   double zg_05 = soft_drop_jet.structure_of<SoftDrop>().symmetry();
   double dr_05 = soft_drop_jet.structure_of<SoftDrop>().delta_R();
   double mu_05 = soft_drop_jet.structure_of<SoftDrop>().mu();
   properties.push_back(MOD::Property("zg_05", zg_05));
   properties.push_back(MOD::Property("dr_05", dr_05));
   properties.push_back(MOD::Property("mu_05", mu_05));

   SoftDrop soft_drop_1(0.0, 0.1);
   PseudoJet soft_drop_jet_1 = soft_drop_1(hardest_jet);
   double zg_1 = soft_drop_jet_1.structure_of<SoftDrop>().symmetry();
   double dr_1 = soft_drop_jet_1.structure_of<SoftDrop>().delta_R();
   double mu_1 = soft_drop_jet_1.structure_of<SoftDrop>().mu();
   properties.push_back(MOD::Property("zg_1", zg_1));  
   properties.push_back(MOD::Property("dr_1", dr_1));  
   properties.push_back(MOD::Property("mu_1", mu_1));  

   SoftDrop soft_drop_2(0.0, 0.2);
   PseudoJet soft_drop_jet_2 = soft_drop_2(hardest_jet);
   double zg_2 = soft_drop_jet_2.structure_of<SoftDrop>().symmetry();
   double dr_2 = soft_drop_jet_2.structure_of<SoftDrop>().delta_R();
   double mu_2 = soft_drop_jet_2.structure_of<SoftDrop>().mu();
   properties.push_back(MOD::Property("zg_2", zg_2));  
   properties.push_back(MOD::Property("dr_2", dr_2));  
   properties.push_back(MOD::Property("mu_2", mu_2));  


   ClusterSequence cs_1(MOD::filter_by_pT(hardest_jet_constituents, 1.00), jet_def_cambridge);

   if (cs_1.inclusive_jets().size() > 0) {
      PseudoJet hardest_jet_pt_1 = cs_1.inclusive_jets()[0];

      SoftDrop soft_drop_pt_1(0.0, 0.05);
      PseudoJet soft_drop_jet_pt_1 = soft_drop_pt_1(hardest_jet_pt_1);
      double zg_05_pt_1 = soft_drop_jet_pt_1.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_05_pt_1", zg_05_pt_1));

      SoftDrop soft_drop_pt_1_1(0.0, 0.1);
      PseudoJet soft_drop_jet_pt_1_1 = soft_drop_pt_1_1(hardest_jet_pt_1);
      double zg_1_pt_1 = soft_drop_jet_pt_1_1.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_1_pt_1", zg_1_pt_1));  

      SoftDrop soft_drop_pt_1_2(0.0, 0.2);
      PseudoJet soft_drop_jet_pt_1_2 = soft_drop_pt_1_2(hardest_jet_pt_1);
      double zg_2_pt_1 = soft_drop_jet_pt_1_2.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_2_pt_1", zg_2_pt_1));  
   }
   else {
      properties.push_back(MOD::Property("zg_05_pt_1", -1.));
      properties.push_back(MOD::Property("zg_1_pt_1", -1.));
      properties.push_back(MOD::Property("zg_2_pt_1", -1.));
   }
   


   
   ClusterSequence cs_2(MOD::filter_by_pT(hardest_jet_constituents, 2.00), jet_def_cambridge);

   if (cs_2.inclusive_jets().size() > 0) {
      PseudoJet hardest_jet_pt_2 = cs_2.inclusive_jets()[0];

      SoftDrop soft_drop_pt_2(0.0, 0.05);
      PseudoJet soft_drop_jet_pt_2 = soft_drop_pt_2(hardest_jet_pt_2);
      double zg_05_pt_2 = soft_drop_jet_pt_2.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_05_pt_2", zg_05_pt_2));

      SoftDrop soft_drop_pt_2_1(0.0, 0.1);
      PseudoJet soft_drop_jet_pt_2_1 = soft_drop_pt_2_1(hardest_jet_pt_2);
      double zg_1_pt_2 = soft_drop_jet_pt_2_1.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_1_pt_2", zg_1_pt_2));  

      SoftDrop soft_drop_pt_2_2(0.0, 0.2);
      PseudoJet soft_drop_jet_pt_2_2 = soft_drop_pt_2_2(hardest_jet_pt_2);
      double zg_2_pt_2 = soft_drop_jet_pt_2_2.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_2_pt_2", zg_2_pt_2));  
   }
   else {
      properties.push_back(MOD::Property("zg_05_pt_2", -1.));
      properties.push_back(MOD::Property("zg_1_pt_2", -1.));
      properties.push_back(MOD::Property("zg_2_pt_2", -1.));
   }
   



   ClusterSequence cs_3(MOD::filter_by_pT(hardest_jet_constituents, 3.00), jet_def_cambridge);

   if (cs_3.inclusive_jets().size() > 0) {
      PseudoJet hardest_jet_pt_3 = cs_3.inclusive_jets()[0];

      SoftDrop soft_drop_pt_3(0.0, 0.05);
      PseudoJet soft_drop_jet_pt_3 = soft_drop_pt_3(hardest_jet_pt_3);
      double zg_05_pt_3 = soft_drop_jet_pt_3.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_05_pt_3", zg_05_pt_3));

      SoftDrop soft_drop_pt_3_1(0.0, 0.1);
      PseudoJet soft_drop_jet_pt_3_1 = soft_drop_pt_3_1(hardest_jet_pt_3);
      double zg_1_pt_3 = soft_drop_jet_pt_3_1.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_1_pt_3", zg_1_pt_3));  

      SoftDrop soft_drop_pt_3_2(0.0, 0.2);
      PseudoJet soft_drop_jet_pt_3_2 = soft_drop_pt_3_2(hardest_jet_pt_3);
      double zg_2_pt_3 = soft_drop_jet_pt_3_2.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_2_pt_3", zg_2_pt_3));  
   }
   else {
      properties.push_back(MOD::Property("zg_05_pt_3", -1.));
      properties.push_back(MOD::Property("zg_1_pt_3", -1.));
      properties.push_back(MOD::Property("zg_2_pt_3", -1.));
   }
   

   

   ClusterSequence cs_5(MOD::filter_by_pT(hardest_jet_constituents, 5.00), jet_def_cambridge);

   if (cs_5.inclusive_jets().size() > 0) {
      PseudoJet hardest_jet_pt_5 = cs_5.inclusive_jets()[0];

      SoftDrop soft_drop_pt_5(0.0, 0.05);
      PseudoJet soft_drop_jet_pt_5 = soft_drop_pt_5(hardest_jet_pt_5);
      double zg_05_pt_5 = soft_drop_jet_pt_5.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_05_pt_5", zg_05_pt_5));

      SoftDrop soft_drop_pt_5_1(0.0, 0.1);
      PseudoJet soft_drop_jet_pt_5_1 = soft_drop_pt_5_1(hardest_jet_pt_5);
      double zg_1_pt_5 = soft_drop_jet_pt_5_1.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_1_pt_5", zg_1_pt_5));  

      SoftDrop soft_drop_pt_5_2(0.0, 0.2);
      PseudoJet soft_drop_jet_pt_5_2 = soft_drop_pt_5_2(hardest_jet_pt_5);
      double zg_2_pt_5 = soft_drop_jet_pt_5_2.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_2_pt_5", zg_2_pt_5));  
   }
   else {
      properties.push_back(MOD::Property("zg_05_pt_5", -1.));
      properties.push_back(MOD::Property("zg_1_pt_5", -1.));
      properties.push_back(MOD::Property("zg_2_pt_5", -1.));
   }
   



   ClusterSequence cs_10(MOD::filter_by_pT(hardest_jet_constituents, 10.00), jet_def_cambridge);

   if (cs_10.inclusive_jets().size() > 0) {
      PseudoJet hardest_jet_pt_10 = cs_10.inclusive_jets()[0];

      SoftDrop soft_drop_pt_10(0.0, 0.05);
      PseudoJet soft_drop_jet_pt_10 = soft_drop_pt_10(hardest_jet_pt_10);
      double zg_05_pt_10 = soft_drop_jet_pt_10.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_05_pt_10", zg_05_pt_10));

      SoftDrop soft_drop_pt_10_1(0.0, 0.1);
      PseudoJet soft_drop_jet_pt_10_1 = soft_drop_pt_10_1(hardest_jet_pt_10);
      double zg_1_pt_10 = soft_drop_jet_pt_10_1.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_1_pt_10", zg_1_pt_10));  

      SoftDrop soft_drop_pt_10_2(0.0, 0.2);
      PseudoJet soft_drop_jet_pt_10_2 = soft_drop_pt_10_2(hardest_jet_pt_10);
      double zg_2_pt_10 = soft_drop_jet_pt_10_2.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_2_pt_10", zg_2_pt_10)); 
   }
   else {
      properties.push_back(MOD::Property("zg_05_pt_10", -1.));
      properties.push_back(MOD::Property("zg_1_pt_10", -1.));
      properties.push_back(MOD::Property("zg_2_pt_10", -1.));
   }

   
   std::vector<fastjet::PseudoJet> charged_constituents = MOD::filter_charged(hardest_jet_constituents);
   
   // Cluster this using Cambridge/Alachen with infinite radius.

   ClusterSequence cs_charged(charged_constituents, jet_def_cambridge);

   if (cs_charged.inclusive_jets().size() > 0 ) {
      PseudoJet hardest_jet_charged = cs_charged.inclusive_jets()[0];

      SoftDrop soft_drop_charged(0.0, 0.05);
      PseudoJet soft_drop_jet_charged = soft_drop_charged(hardest_jet_charged);
      double zg_charged_05 = soft_drop_jet_charged.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_charged_05", zg_charged_05));

      SoftDrop soft_drop_charged_2(0.0, 0.1);
      PseudoJet soft_drop_jet_charged_2 = soft_drop_charged_2(hardest_jet_charged);
      double zg_charged_1 = soft_drop_jet_charged_2.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_charged_1", zg_charged_1));  

      SoftDrop soft_drop_charged_3(0.0, 0.2);
      PseudoJet soft_drop_jet_charged_3 = soft_drop_charged_3(hardest_jet_charged);
      double zg_charged_2 = soft_drop_jet_charged_3.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_charged_2", zg_charged_2));     
   }
   else {
      properties.push_back(MOD::Property("zg_charged_05", -1.0));
      properties.push_back(MOD::Property("zg_charged_05", -1.0));
      properties.push_back(MOD::Property("zg_charged_05", -1.0));      
   }
   

   // Stuff before and after SoftDrop.
   

   // Hardest Jet pT before and after SoftDrop.
   properties.push_back(MOD::Property("pT_after_SD", soft_drop_jet.pt()));

   // Jet mass and constituent multiplicity before and after SoftDrop.

   properties.push_back( MOD::Property("mul_pre_SD", (int) hardest_jet_constituents.size()) );
   properties.push_back( MOD::Property("mul_post_SD", (int) soft_drop(hardest_jet).constituents().size() ) );   

   properties.push_back( MOD::Property("mass_pre_SD", hardest_jet.m()) );
   properties.push_back( MOD::Property("mass_post_SD", soft_drop(hardest_jet).m()) );

   
   // Jet mass and multiplicity before and after SD for charged particles only.

   if (cs_charged.inclusive_jets().size() > 0 ) {
      PseudoJet hardest_jet_charged = cs_charged.inclusive_jets()[0];

      properties.push_back( MOD::Property("chrg_mul_pre_SD", (int) hardest_jet_charged.constituents().size()) );
      properties.push_back( MOD::Property("chrg_mul_post_SD", (int) soft_drop(hardest_jet_charged).constituents().size()) );

      properties.push_back( MOD::Property("chrg_mass_pre_SD", hardest_jet_charged.m()) );
      properties.push_back( MOD::Property("chrg_mass_post_SD", soft_drop(hardest_jet_charged).m()) );

      SoftDrop chrg_soft_drop_05(0.0, 0.05);
      SoftDrop chrg_soft_drop_1(0.0, 0.1);
      SoftDrop chrg_soft_drop_2(0.0, 0.2);

      properties.push_back( MOD::Property("chrg_dr_05", chrg_soft_drop_05(hardest_jet_charged).structure_of<SoftDrop>().delta_R()));
      properties.push_back( MOD::Property("chrg_dr_1", chrg_soft_drop_1(hardest_jet_charged).structure_of<SoftDrop>().delta_R()));
      properties.push_back( MOD::Property("chrg_dr_2", chrg_soft_drop_2(hardest_jet_charged).structure_of<SoftDrop>().delta_R()));

   }
   else {
      properties.push_back( MOD::Property("chrg_mul_pre_SD", -1. ));
      properties.push_back( MOD::Property("chrg_mul_post_SD", -1. ));

      properties.push_back( MOD::Property("chrg_mass_pre_SD", -1. ));
      properties.push_back( MOD::Property("chrg_mass_post_SD", -1. ));

      properties.push_back( MOD::Property("chrg_dr_05", -1. ));
      properties.push_back( MOD::Property("chrg_dr_1", -1. ));
      properties.push_back( MOD::Property("chrg_dr_2", -1. ));
   }
   


   properties.push_back( MOD::Property("fra_energy_loss", (hardest_jet.E() - soft_drop(hardest_jet).E()) / hardest_jet.E() ) );
   properties.push_back( MOD::Property("hardest_eta", hardest_jet.eta()) );
   
   


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