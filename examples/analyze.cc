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

      // Write out version info in the output file for the "syncing plots" thing to work (as it needs to figure out which directory to put things into).
      if (event_serial_number == 1)
         output_file << "%" << " Version " << event_being_read.version() << endl;

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


   MOD::CalibratedJet hardest_uncorrected_jet = event_being_read.hardest_jet(true, false, true, "loose", 2.4);
   MOD::CalibratedJet hardest_corrected_jet = event_being_read.hardest_jet(true, true, true, "loose", 2.4);

   if (event_being_read.CMS_jets().size() == 0) {
      return;
   }
   

   vector<MOD::Property> properties;

   properties.push_back(MOD::Property("# Entry", "  Entry"));

   properties.push_back(MOD::Property("Event_Number", event_being_read.event_number()));
   properties.push_back(MOD::Property("Run_Number", event_being_read.run_number()));
   


   properties.push_back(MOD::Property("Uncor_Hardest_pT", hardest_uncorrected_jet.pseudojet().pt()));
   properties.push_back(MOD::Property("Cor_Hardest_pT", hardest_corrected_jet.pseudojet().pt()));

   properties.push_back(MOD::Property("Prescale", event_being_read.assigned_trigger_prescale()));
   properties.push_back(MOD::Property("Trigger_Name", event_being_read.assigned_trigger_name()));


   // Run AK5 clustering with FastJet to get zg value.

   // JetDefinition jet_def(antikt_algorithm, 0.5);
   // ClusterSequence cs(event_being_read.pseudojets(), jet_def);
   // vector<PseudoJet> ak5_jets = sorted_by_pt(cs.inclusive_jets(3.0));


   // PseudoJet hardest_jet = ak5_jets[0];
   // // hardest_jet *= event_being_read.hardest_jet_JEC();
   PseudoJet hardest_jet_fo

   double beta = 0;

   SoftDrop soft_drop(beta, 0.05);
   PseudoJet soft_drop_jet = soft_drop(hardest_jet);
   double zg_05 = soft_drop_jet.structure_of<SoftDrop>().symmetry();
   double dr_05 = soft_drop_jet.structure_of<SoftDrop>().delta_R();
   double mu_05 = soft_drop_jet.structure_of<SoftDrop>().mu();
   properties.push_back(MOD::Property("zg_05", zg_05));
   properties.push_back(MOD::Property("dr_05", dr_05));
   properties.push_back(MOD::Property("mu_05", mu_05));

   SoftDrop soft_drop_1(beta, 0.1);
   PseudoJet soft_drop_jet_1 = soft_drop_1(hardest_jet);
   double zg_1 = soft_drop_jet_1.structure_of<SoftDrop>().symmetry();
   double dr_1 = soft_drop_jet_1.structure_of<SoftDrop>().delta_R();
   double mu_1 = soft_drop_jet_1.structure_of<SoftDrop>().mu();
   properties.push_back(MOD::Property("zg_1", zg_1));  
   properties.push_back(MOD::Property("dr_1", dr_1));  
   properties.push_back(MOD::Property("mu_1", mu_1));  

   SoftDrop soft_drop_2(beta, 0.2);
   PseudoJet soft_drop_jet_2 = soft_drop_2(hardest_jet);
   double zg_2 = soft_drop_jet_2.structure_of<SoftDrop>().symmetry();
   double dr_2 = soft_drop_jet_2.structure_of<SoftDrop>().delta_R();
   double mu_2 = soft_drop_jet_2.structure_of<SoftDrop>().mu();
   properties.push_back(MOD::Property("zg_2", zg_2));  
   properties.push_back(MOD::Property("dr_2", dr_2));  
   properties.push_back(MOD::Property("mu_2", mu_2));  


   properties.push_back(MOD::Property("Hardest_PFC_pdgId", event_being_read.hardest_pfcandidate().pdgId()));
   properties.push_back(MOD::Property("Hardest_PFC_pT", event_being_read.hardest_pfcandidate().pseudojet().pt()));




   // Run AK5 clustering with FastJet to get zg value.

   JetDefinition jet_def_pt_1(antikt_algorithm, 0.5);
   ClusterSequence cs_pt_1(event_being_read.pseudojets(1.00), jet_def_pt_1);
   vector<PseudoJet> ak5_jets_pt_1 = sorted_by_pt(cs_pt_1.inclusive_jets(3.0));


   PseudoJet hardest_jet_pt_1 = ak5_jets_pt_1[0];

   double beta_pt_1 = 0;

   SoftDrop soft_drop_pt_1(beta_pt_1, 0.05);
   PseudoJet soft_drop_jet_pt_1 = soft_drop_pt_1(hardest_jet_pt_1);
   double zg_05_pt_1 = soft_drop_jet_pt_1.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_05_pt_1", zg_05_pt_1));

   SoftDrop soft_drop_pt_1_1(beta_pt_1, 0.1);
   PseudoJet soft_drop_jet_pt_1_1 = soft_drop_pt_1_1(hardest_jet_pt_1);
   double zg_1_pt_1 = soft_drop_jet_pt_1_1.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_1_pt_1", zg_1_pt_1));  

   SoftDrop soft_drop_pt_1_2(beta_pt_1, 0.2);
   PseudoJet soft_drop_jet_pt_1_2 = soft_drop_pt_1_2(hardest_jet_pt_1);
   double zg_2_pt_1 = soft_drop_jet_pt_1_2.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_2_pt_1", zg_2_pt_1));  



   // Run AK5 clustering with FastJet to get zg value.

   JetDefinition jet_def_pt_2(antikt_algorithm, 0.5);
   ClusterSequence cs_pt_2(event_being_read.pseudojets(2.00), jet_def_pt_2);
   vector<PseudoJet> ak5_jets_pt_2 = sorted_by_pt(cs_pt_2.inclusive_jets(3.0));

   PseudoJet hardest_jet_pt_2 = ak5_jets_pt_2[0];

   double beta_pt_2 = 0;

   SoftDrop soft_drop_pt_2(beta_pt_2, 0.05);
   PseudoJet soft_drop_jet_pt_2 = soft_drop_pt_2(hardest_jet_pt_2);
   double zg_05_pt_2 = soft_drop_jet_pt_2.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_05_pt_2", zg_05_pt_2));

   SoftDrop soft_drop_pt_2_1(beta_pt_2, 0.1);
   PseudoJet soft_drop_jet_pt_2_1 = soft_drop_pt_2_1(hardest_jet_pt_2);
   double zg_1_pt_2 = soft_drop_jet_pt_2_1.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_1_pt_2", zg_1_pt_2));  

   SoftDrop soft_drop_pt_2_2(beta_pt_2, 0.2);
   PseudoJet soft_drop_jet_pt_2_2 = soft_drop_pt_2_2(hardest_jet_pt_2);
   double zg_2_pt_2 = soft_drop_jet_pt_2_2.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_2_pt_2", zg_2_pt_2));  


   // Run AK5 clustering with FastJet to get zg value.

   JetDefinition jet_def_pt_3(antikt_algorithm, 0.5);
   ClusterSequence cs_pt_3(event_being_read.pseudojets(3.00), jet_def_pt_3);
   vector<PseudoJet> ak5_jets_pt_3 = sorted_by_pt(cs_pt_3.inclusive_jets(3.0));

   PseudoJet hardest_jet_pt_3 = ak5_jets_pt_3[0];

   double beta_pt_3 = 0;

   SoftDrop soft_drop_pt_3(beta_pt_3, 0.05);
   PseudoJet soft_drop_jet_pt_3 = soft_drop_pt_3(hardest_jet_pt_3);
   double zg_05_pt_3 = soft_drop_jet_pt_3.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_05_pt_3", zg_05_pt_3));

   SoftDrop soft_drop_pt_3_1(beta_pt_3, 0.1);
   PseudoJet soft_drop_jet_pt_3_1 = soft_drop_pt_3_1(hardest_jet_pt_3);
   double zg_1_pt_3 = soft_drop_jet_pt_3_1.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_1_pt_3", zg_1_pt_3));  

   SoftDrop soft_drop_pt_3_2(beta_pt_3, 0.2);
   PseudoJet soft_drop_jet_pt_3_2 = soft_drop_pt_3_2(hardest_jet_pt_3);
   double zg_2_pt_3 = soft_drop_jet_pt_3_2.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_2_pt_3", zg_2_pt_3));  

   // Run AK5 clustering with FastJet to get zg value.

   JetDefinition jet_def_pt_5(antikt_algorithm, 0.5);
   ClusterSequence cs_pt_5(event_being_read.pseudojets(5.00), jet_def_pt_5);
   vector<PseudoJet> ak5_jets_pt_5 = sorted_by_pt(cs_pt_5.inclusive_jets(3.0));


   PseudoJet hardest_jet_pt_5 = ak5_jets_pt_5[0];

   double beta_pt_5 = 0;

   SoftDrop soft_drop_pt_5(beta_pt_5, 0.05);
   PseudoJet soft_drop_jet_pt_5 = soft_drop_pt_5(hardest_jet_pt_5);
   double zg_05_pt_5 = soft_drop_jet_pt_5.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_05_pt_5", zg_05_pt_5));

   SoftDrop soft_drop_pt_5_1(beta_pt_5, 0.1);
   PseudoJet soft_drop_jet_pt_5_1 = soft_drop_pt_5_1(hardest_jet_pt_5);
   double zg_1_pt_5 = soft_drop_jet_pt_5_1.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_1_pt_5", zg_1_pt_5));  

   SoftDrop soft_drop_pt_5_2(beta_pt_5, 0.2);
   PseudoJet soft_drop_jet_pt_5_2 = soft_drop_pt_5_2(hardest_jet_pt_5);
   double zg_2_pt_5 = soft_drop_jet_pt_5_2.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_2_pt_5", zg_2_pt_5));  


   // Run AK5 clustering with FastJet to get zg value.

   JetDefinition jet_def_pt_10(antikt_algorithm, 0.5);
   ClusterSequence cs_pt_10(event_being_read.pseudojets(10.00), jet_def_pt_10);
   vector<PseudoJet> ak5_jets_pt_10 = sorted_by_pt(cs_pt_10.inclusive_jets(3.0));


   PseudoJet hardest_jet_pt_10 = ak5_jets_pt_10[0];

   double beta_pt_10 = 0;

   SoftDrop soft_drop_pt_10(beta_pt_10, 0.05);
   PseudoJet soft_drop_jet_pt_10 = soft_drop_pt_10(hardest_jet_pt_10);
   double zg_05_pt_10 = soft_drop_jet_pt_10.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_05_pt_10", zg_05_pt_10));

   SoftDrop soft_drop_pt_10_1(beta_pt_10, 0.1);
   PseudoJet soft_drop_jet_pt_10_1 = soft_drop_pt_10_1(hardest_jet_pt_10);
   double zg_1_pt_10 = soft_drop_jet_pt_10_1.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_1_pt_10", zg_1_pt_10));  

   SoftDrop soft_drop_pt_10_2(beta_pt_10, 0.2);
   PseudoJet soft_drop_jet_pt_10_2 = soft_drop_pt_10_2(hardest_jet_pt_10);
   double zg_2_pt_10 = soft_drop_jet_pt_10_2.structure_of<SoftDrop>().symmetry();
   properties.push_back(MOD::Property("zg_2_pt_10", zg_2_pt_10));  





   vector<PseudoJet> hardest_jet_constituents = ak5_jets[0].constituents();

   vector<PseudoJet> hardest_jet_charged_constituents;

   for (unsigned i = 0; i < hardest_jet_constituents.size(); i++) {
      if ( (abs(hardest_jet_constituents[i].user_index()) == 211) || (abs(hardest_jet_constituents[i].user_index()) == 11) || (abs(hardest_jet_constituents[i].user_index()) == 13) ) {
         hardest_jet_charged_constituents.push_back(hardest_jet_constituents[i]);
      }
   }

   JetDefinition jet_def_charged(aachen_algorithm, 1000);
   ClusterSequence cs_charged(hardest_jet_charged_constituents, jet_def_charged);
   vector<PseudoJet> ak5_jet_with_charged_particles_only = sorted_by_pt(cs_charged.inclusive_jets(3.0));

   if (ak5_jet_with_charged_particles_only.size() > 0) {
      PseudoJet hardest_jet_charged = ak5_jet_with_charged_particles_only[0];

      double beta_charged = 0;

      SoftDrop soft_drop_charged(beta_charged, 0.05);
      PseudoJet soft_drop_jet_charged = soft_drop_charged(hardest_jet_charged);
      double zg_charged_05 = soft_drop_jet_charged.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_charged_05", zg_charged_05));

      SoftDrop soft_drop_charged_2(beta_charged, 0.1);
      PseudoJet soft_drop_jet_charged_2 = soft_drop_charged_2(hardest_jet_charged);
      double zg_charged_1 = soft_drop_jet_charged_2.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_charged_1", zg_charged_1));  

      SoftDrop soft_drop_charged_3(beta_charged, 0.2);
      PseudoJet soft_drop_jet_charged_3 = soft_drop_charged_3(hardest_jet_charged);
      double zg_charged_2 = soft_drop_jet_charged_3.structure_of<SoftDrop>().symmetry();
      properties.push_back(MOD::Property("zg_charged_2", zg_charged_2));  
   }
   else {
      properties.push_back(MOD::Property("zg_charged_05", -1.00));      
      properties.push_back(MOD::Property("zg_charged_1", -1.00));      
      properties.push_back(MOD::Property("zg_charged_2", -1.00));  
   }




   // Stuff before and after SoftDrop.
   

   // Hardest Jet pT before and after SoftDrop.
   properties.push_back(MOD::Property("pT_after_SD", soft_drop_jet.pt()));

   // Constituent multiplicity before and after SoftDrop.

   properties.push_back( MOD::Property("multiplicity_before_SD", (int) hardest_jet.constituents().size()) );
   properties.push_back( MOD::Property("multiplicity_after_SD", (int) soft_drop(hardest_jet).constituents().size()) );
   

   // Must be "hardest_uncorrected_jet" because the corrected jet has a JEC set to 1.
   properties.push_back( MOD::Property("jec", hardest_uncorrected_jet.JEC()) );
   
   properties.push_back( MOD::Property("jet_area", hardest_uncorrected_jet.area()) );

   

   properties.push_back( MOD::Property("jet_mass_before_SD", hardest_jet.m()) );
   properties.push_back( MOD::Property("jet_mass_after_SD", soft_drop(hardest_jet).m()) );

   properties.push_back( MOD::Property("fractional_energy_loss", (hardest_jet.E() - soft_drop(hardest_jet).E()) / hardest_jet.E() ) );
   
   properties.push_back( MOD::Property("hardest_eta", hardest_jet.eta()) );


   // Jet mass and multiplicity before and after SD for charged particles only.
   PseudoJet hardest_jet_charged = ak5_jet_with_charged_particles_only[0];

   properties.push_back( MOD::Property("charged_multi_before_SD", (int) hardest_jet_charged.constituents().size()) );
   properties.push_back( MOD::Property("charged_multi_after_SD", (int) soft_drop(hardest_jet_charged).constituents().size()) );

   properties.push_back( MOD::Property("charged_jet_mass_before_SD", hardest_jet_charged.m()) );
   properties.push_back( MOD::Property("charged_jet_mass_after_SD", soft_drop(hardest_jet_charged).m()) );

   

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