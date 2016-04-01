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
#include "../interface/property.h"

using namespace std;
using namespace fastjet;
using namespace contrib;

void analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number);

double angularity_lambda(PseudoJet jet, int k, int beta);

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
      
      if( (event_serial_number % 1000) == 0 )
         cout << "Processing event number " << event_serial_number << endl;

      // Write out version info in the output file for the "syncing plots" thing to work (as it needs to figure out which directory to put things into).
      if (event_serial_number == 1)
         output_file << "%" << " Version " << event_being_read.version() << endl;

      
      analyze_event(event_being_read, output_file, event_serial_number);
      
      event_being_read = MOD::Event();
      event_serial_number++;

   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_event(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number) {

   JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);

   PseudoJet hardest_jet = event_being_read.hardest_jet();
   vector<PseudoJet> hardest_jet_constituents = hardest_jet.constituents();
   
   SoftDrop soft_drop(0.0, 0.1);
   PseudoJet soft_drop_jet = soft_drop(hardest_jet);

   
   vector<MOD::Property> properties;
   
   properties.push_back(MOD::Property("# Entry", "  Entry"));

   properties.push_back(MOD::Property("prescale", event_being_read.weight()));
   properties.push_back(MOD::Property("hardest_pT", hardest_jet.pt()));

   properties.push_back( MOD::Property("frac_pT_loss", (hardest_jet.pt() - soft_drop(hardest_jet).pt()) / hardest_jet.pt() ) );
   properties.push_back( MOD::Property("hardest_eta", hardest_jet.eta()) );
   properties.push_back( MOD::Property("hardest_phi", hardest_jet.phi()) );


   vector<pair<string, double>> zg_cuts { make_pair("05", 0.05), make_pair("10", 0.1), make_pair("20", 0.2) };

   // zg, dr, and mu for zg_cuts of 0.05, 0.1 and 0.2.

   for (unsigned i = 0; i < zg_cuts.size(); i++) {

      string label = zg_cuts[i].first;
      double zg_cut = zg_cuts[i].second;

      SoftDrop soft_drop(0.0, zg_cut);
   
      PseudoJet soft_drop_jet = soft_drop(hardest_jet);

      properties.push_back(MOD::Property("zg_" + label, soft_drop_jet.structure_of<SoftDrop>().symmetry()));
      properties.push_back(MOD::Property("Rg_" + label, soft_drop_jet.structure_of<SoftDrop>().delta_R()));
      properties.push_back(MOD::Property("mu_" + label, soft_drop_jet.structure_of<SoftDrop>().mu()));
   }

   // SoftKiller.
   // pT cut of 1, 2, 3, 5, 10 GeV used for SoftKiller.

   vector<int> pT_cuts {1, 2, 3, 5, 10};

   for (unsigned i = 0; i < pT_cuts.size(); i++) {
      int pT_cut = pT_cuts[i];

      Selector pT_selector = SelectorPtMin(pT_cut);
      ClusterSequence cs = ClusterSequence(pT_selector(hardest_jet_constituents), jet_def_cambridge);

      if (cs.inclusive_jets().size() > 0) {
         
         PseudoJet hardest_jet_pT_cut = cs.inclusive_jets()[0];

         for (unsigned j = 0; j < zg_cuts.size(); j++) {
            string label = zg_cuts[j].first;
            double zg_cut = zg_cuts[j].second;

            SoftDrop soft_drop_pT(0.0, zg_cut);
            PseudoJet soft_drop_jet_pT = soft_drop_pT(hardest_jet_pT_cut);
            
            properties.push_back(MOD::Property("zg_" + label + "_pT_" + to_string(pT_cut), soft_drop_jet_pT.structure_of<SoftDrop>().symmetry()));
            properties.push_back(MOD::Property("Rg_" + label + "_pT_" + to_string(pT_cut), soft_drop_jet_pT.structure_of<SoftDrop>().delta_R()));
            properties.push_back(MOD::Property("mu_" + label + "_pT_" + to_string(pT_cut), soft_drop_jet_pT.structure_of<SoftDrop>().mu()));

         }
      }
      else {
         for (unsigned j = 0; j < zg_cuts.size(); j++) {
            string label = zg_cuts[j].first;
            properties.push_back(MOD::Property("zg_" + label + "_pT_" + to_string(pT_cut), -1.));
            properties.push_back(MOD::Property("Rg_" + label + "_pT_" + to_string(pT_cut), -1.));
            properties.push_back(MOD::Property("mu_" + label + "_pT_" + to_string(pT_cut), -1.));
         }
      }
   }   

   // Analysis related to the effects of SoftDrop- observables before and after SoftDrop.
   
   properties.push_back(MOD::Property("pT_after_SD", soft_drop_jet.pt()));

   properties.push_back( MOD::Property("mul_pre_SD", (int) hardest_jet_constituents.size()) );
   properties.push_back( MOD::Property("mul_post_SD", (int) soft_drop(hardest_jet).constituents().size() ) );   

   properties.push_back( MOD::Property("mass_pre_SD", hardest_jet.m()) );
   properties.push_back( MOD::Property("mass_post_SD", soft_drop(hardest_jet).m()) );


   properties.push_back( MOD::Property("pT_D_pre_SD", angularity_lambda(hardest_jet, 2, 0)) );
   properties.push_back( MOD::Property("pT_D_post_SD", angularity_lambda(soft_drop(hardest_jet), 2, 0)) );

   properties.push_back( MOD::Property("LHA_pre_SD", angularity_lambda(hardest_jet, 1, 0.5)) );
   properties.push_back( MOD::Property("LHA_post_SD", angularity_lambda(soft_drop(hardest_jet), 1, 0.5)) );

   properties.push_back( MOD::Property("width_pre_SD", angularity_lambda(hardest_jet, 1, 1)) );
   properties.push_back( MOD::Property("width_post_SD", angularity_lambda(soft_drop(hardest_jet), 1, 1)) );

   properties.push_back( MOD::Property("thrust_pre_SD", angularity_lambda(hardest_jet, 1, 2)) );
   properties.push_back( MOD::Property("thrust_post_SD", angularity_lambda(soft_drop(hardest_jet), 1, 2)) );






   





   // ================================================================ Track Based Analysis ================================================================





   Selector pT_1_GeV_selector = SelectorPtMin(1.0);

   // Get all charged particles with 1 GeV particles removed.
   std::vector<fastjet::PseudoJet> charged_1GeV_constituents = pT_1_GeV_selector(MOD::filter_charged(hardest_jet_constituents));

   // Cluster them using Cambridge/Alachen with infinite radius. This makes sure that we get the same jets as "regular" ak5 jets except now with just charged particles.
   ClusterSequence cs_charged_1_GeV(charged_1GeV_constituents, jet_def_cambridge);

   if (cs_charged_1_GeV.inclusive_jets().size() > 0 ) {

      PseudoJet hardest_charged_1_GeV_removed_jet = cs_charged_1_GeV.inclusive_jets()[0];

      for (unsigned i = 0; i < zg_cuts.size(); i++) {
         string label = zg_cuts[i].first;
         double zg_cut = zg_cuts[i].second;

         SoftDrop soft_drop_charged(0.0, zg_cut);
         PseudoJet soft_drop_jet_charged = soft_drop_charged(hardest_charged_1_GeV_removed_jet);

         properties.push_back(MOD::Property("charged_1_GeV_zg_" + label, soft_drop_jet_charged.structure_of<SoftDrop>().symmetry()));
         properties.push_back(MOD::Property("charged_1_GeV_Rg_" + label, soft_drop_jet_charged.structure_of<SoftDrop>().delta_R()));
         properties.push_back(MOD::Property("charged_1_GeV_mu_" + label, soft_drop_jet_charged.structure_of<SoftDrop>().mu()));
      }
   }
   else {
      for (unsigned i = 0; i < zg_cuts.size(); i++) {
         string label = zg_cuts[i].first;
         properties.push_back(MOD::Property("charged_1_GeV_zg_" + label, -1.0));
         properties.push_back(MOD::Property("charged_1_GeV_Rg_" + label, -1.0));
         properties.push_back(MOD::Property("charged_1_GeV_mu_" + label, -1.0));
      }
   }


   // Before and after SoftDrop for 1 GeV removed charged particles only.

   if (cs_charged_1_GeV.inclusive_jets().size() > 0 ) {
      PseudoJet hardest_charged_1_GeV_removed_jet = cs_charged_1_GeV.inclusive_jets()[0];

      properties.push_back( MOD::Property("chrg_1_GeV_mul_pre_SD", (int) hardest_charged_1_GeV_removed_jet.constituents().size()) );
      properties.push_back( MOD::Property("chrg_1_GeV_mul_post_SD", (int) soft_drop(hardest_charged_1_GeV_removed_jet).constituents().size()) );

      properties.push_back( MOD::Property("chrg_1_GeV_mass_pre_SD", hardest_charged_1_GeV_removed_jet.m()) );
      properties.push_back( MOD::Property("chrg_1_GeV_mass_post_SD", soft_drop(hardest_charged_1_GeV_removed_jet).m()) );


      properties.push_back( MOD::Property("chrg_1_GeV_pT_D_pre_SD", angularity_lambda(hardest_charged_1_GeV_removed_jet, 2, 0)) );
      properties.push_back( MOD::Property("chrg_1_GeV_pT_D_post_SD", angularity_lambda(soft_drop(hardest_charged_1_GeV_removed_jet), 2, 0)) );

      properties.push_back( MOD::Property("chrg_1_GeV_LHA_pre_SD", angularity_lambda(hardest_charged_1_GeV_removed_jet, 1, 0.5)) );
      properties.push_back( MOD::Property("chrg_1_GeV_LHA_post_SD", angularity_lambda(soft_drop(hardest_charged_1_GeV_removed_jet), 1, 0.5)) );

      properties.push_back( MOD::Property("chrg_1_GeV_width_pre_SD", angularity_lambda(hardest_charged_1_GeV_removed_jet, 1, 1)) );
      properties.push_back( MOD::Property("chrg_1_GeV_width_post_SD", angularity_lambda(soft_drop(hardest_charged_1_GeV_removed_jet), 1, 1)) );

      properties.push_back( MOD::Property("chrg_1_GeV_thrust_pre_SD", angularity_lambda(hardest_charged_1_GeV_removed_jet, 1, 2)) );
      properties.push_back( MOD::Property("chrg_1_GeV_thrust_post_SD", angularity_lambda(soft_drop(hardest_charged_1_GeV_removed_jet), 1, 2)) );



   }
   else {
      properties.push_back( MOD::Property("chrg_1_GeV_mul_pre_SD", -1. ));
      properties.push_back( MOD::Property("chrg_1_GeV_mul_post_SD", -1. ));

      properties.push_back( MOD::Property("chrg_1_GeV_mass_pre_SD", -1. ));
      properties.push_back( MOD::Property("chrg_1_GeV_mass_post_SD", -1. ));


      properties.push_back( MOD::Property("chrg_1_GeV_pT_D_pre_SD", -1.) );
      properties.push_back( MOD::Property("chrg_1_GeV_pT_D_post_SD", -1.) );

      properties.push_back( MOD::Property("chrg_1_GeV_LHA_pre_SD", -1.) );
      properties.push_back( MOD::Property("chrg_1_GeV_LHA_post_SD", -1.) );

      properties.push_back( MOD::Property("chrg_1_GeV_width_pre_SD", -1.) );
      properties.push_back( MOD::Property("chrg_1_GeV_width_post_SD", -1.) );

      properties.push_back( MOD::Property("chrg_1_GeV_thrust_pre_SD", -1.) );
      properties.push_back( MOD::Property("chrg_1_GeV_thrust_post_SD", -1.) );
   }



   




   

   
   



   

   // Now that we've calculated all observables, write them out.

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





double angularity_lambda(PseudoJet jet, int k, int beta) {
   
   double lambda = 0.0;

   double R = 0.5;   // Jet Radius.

   double total_pT = 0.0;
   for (unsigned j = 0; j < jet.constituents().size(); j++) {
      total_pT += jet.constituents()[j].pt();
   }

   for (unsigned i = 0; i < jet.constituents().size(); i++) {
      
      PseudoJet constituent = jet.constituents()[i];

      double z_i = constituent.pt() / total_pT;
      
      double delta_R = pow( pow(constituent.rap(), 2) + pow(constituent.phi(), 2), 0.5 );

      double theta_i = delta_R / R;

      lambda += pow(z_i, k) * pow(theta_i, beta);

   }

   return lambda;

}