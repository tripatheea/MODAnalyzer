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

void analyze_pfc(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number);


void write_output(int & event_serial_number, ofstream & output_file, PseudoJet vector_1, PseudoJet vector_2, float weight);

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

      analyze_pfc(event_being_read, output_file, event_serial_number);
      
      event_being_read = MOD::Event();
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_pfc(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number) {

   float weight = event_being_read.weight();
  
   JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);

   PseudoJet hardest_jet = event_being_read.hardest_jet();

   vector<PseudoJet> hardest_jet_pfcs = hardest_jet.constituents();

   // electron, muon, pion, photon.
   vector<int> electrons = {-11, 11};
   vector<int> muons = {-13, 13};
   vector<int> pions = {-211, 211};
   vector<int> photons = {22, 22};

   for (unsigned i = 0; i < hardest_jet_pfcs.size(); i++) {

		// vector<MOD::Property> properties;


      // Electron.

      // - - and - +.

      if ( hardest_jet_pfcs[i].user_info<MOD::InfoPFC>().pdgId() == electrons[0] ) {
         
         for (unsigned j = 0; j < hardest_jet_pfcs.size(); j++) {

            // - +
            
            if ( (i != j) & (hardest_jet_pfcs[j].user_info<MOD::InfoPFC>().pdgId() == electrons[1])) {
               write_output(event_serial_number, output_file, hardest_jet_pfcs[i], hardest_jet_pfcs[j], weight);
               event_serial_number++;
            }
         }
      } 


      // + - and + +.

      if ( hardest_jet_pfcs[i].user_info<MOD::InfoPFC>().pdgId() == electrons[1] ) {
         
         for (unsigned j = 0; j < hardest_jet_pfcs.size(); j++) {

            // + -
            
            if ( (i != j) & (hardest_jet_pfcs[j].user_info<MOD::InfoPFC>().pdgId() == electrons[0])) {
               write_output(event_serial_number, output_file, hardest_jet_pfcs[i], hardest_jet_pfcs[j], weight);
               event_serial_number++;
            }
         }
      }



      // Muon.

      // - - and - +.

      if ( hardest_jet_pfcs[i].user_info<MOD::InfoPFC>().pdgId() == muons[0] ) {
         
         for (unsigned j = 0; j < hardest_jet_pfcs.size(); j++) {

            // - +
            
            if ( (i != j) & (hardest_jet_pfcs[j].user_info<MOD::InfoPFC>().pdgId() == muons[1])) {
               write_output(event_serial_number, output_file, hardest_jet_pfcs[i], hardest_jet_pfcs[j], weight);
               event_serial_number++;
            }
         }
      } 


      // + - and + +.

      if ( hardest_jet_pfcs[i].user_info<MOD::InfoPFC>().pdgId() == muons[1] ) {
         
         for (unsigned j = 0; j < hardest_jet_pfcs.size(); j++) {

            // + -
            
            if ( (i != j) & (hardest_jet_pfcs[j].user_info<MOD::InfoPFC>().pdgId() == muons[0])) {
               write_output(event_serial_number, output_file, hardest_jet_pfcs[i], hardest_jet_pfcs[j], weight);
               event_serial_number++;
            }
         }
      }





      // Pion.
      
      // - - and - +.

      if ( hardest_jet_pfcs[i].user_info<MOD::InfoPFC>().pdgId() == pions[0] ) {
         
         for (unsigned j = 0; j < hardest_jet_pfcs.size(); j++) {

            if ( (i != j) & (hardest_jet_pfcs[j].user_info<MOD::InfoPFC>().pdgId() == pions[1])) {          // - +
               write_output(event_serial_number, output_file, hardest_jet_pfcs[i], hardest_jet_pfcs[j], weight);
               event_serial_number++;
            }
            else if ( (i != j) & (hardest_jet_pfcs[j].user_info<MOD::InfoPFC>().pdgId() == pions[0])) {     // - - 
               write_output(event_serial_number, output_file, hardest_jet_pfcs[i], hardest_jet_pfcs[j], weight);
               event_serial_number++;
            }
         
         }

      } 


      // + - and + +.

      if ( hardest_jet_pfcs[i].user_info<MOD::InfoPFC>().pdgId() == pions[1] ) {
         
         for (unsigned j = 0; j < hardest_jet_pfcs.size(); j++) {

            if ( (i != j) & (hardest_jet_pfcs[j].user_info<MOD::InfoPFC>().pdgId() == pions[0])) {          // + -
               write_output(event_serial_number, output_file, hardest_jet_pfcs[i], hardest_jet_pfcs[j], weight);
               event_serial_number++;
            }
            else if ( (i != j) & (hardest_jet_pfcs[j].user_info<MOD::InfoPFC>().pdgId() == pions[1])) {          // + +
               write_output(event_serial_number, output_file, hardest_jet_pfcs[i], hardest_jet_pfcs[j], weight);
               event_serial_number++;
            }
         }
      }








      // Photon.

      if ( hardest_jet_pfcs[i].user_info<MOD::InfoPFC>().pdgId() == photons[0] ) {
         
         for (unsigned j = 0; j < hardest_jet_pfcs.size(); j++) {
            
            if ( (i != j) & (hardest_jet_pfcs[j].user_info<MOD::InfoPFC>().pdgId() == photons[1])) {
               write_output(event_serial_number, output_file, hardest_jet_pfcs[i], hardest_jet_pfcs[j], weight);
               event_serial_number++;
            }
         }
      } 

      
   }
}



void write_output(int & event_serial_number, ofstream & output_file, PseudoJet vector_1, PseudoJet vector_2, float weight) {
   

   vector<MOD::Property> properties;

   properties.push_back(MOD::Property("# Entry", "  Entry"));
   properties.push_back(MOD::Property("prescale", weight));
   properties.push_back(MOD::Property("pdg_id_1", vector_1.user_info<MOD::InfoPFC>().pdgId() ));
   properties.push_back(MOD::Property("pdg_id_2", vector_2.user_info<MOD::InfoPFC>().pdgId() ));
   properties.push_back(MOD::Property("invariant_mass", (vector_1 + vector_2).m() ));



   string name;
   
   int padding = 15;


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

