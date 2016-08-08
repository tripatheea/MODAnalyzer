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


#include "../interface/Event.h"
#include "../interface/Property.h"

using namespace std;
using namespace fastjet;
using namespace contrib;

void analyze_pfc(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number);

double angularity_lambda(PseudoJet jet, float k, float beta);

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

      
      analyze_pfc(event_being_read, output_file, event_serial_number);
      
      event_being_read = MOD::Event();
      event_serial_number++;

   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_pfc(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number) {

  
   JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);

   PseudoJet hardest_jet = event_being_read.hardest_jet();

   vector<PseudoJet> hardest_jet_pfcs = hardest_jet.constituents();

   for (unsigned i = 0; i < hardest_jet_pfcs.size(); i++) {

   		vector<MOD::Property> properties;

   		properties.push_back(MOD::Property("# Entry", "  Entry"));

         properties.push_back(MOD::Property("prescale", event_being_read.weight()));
	   	properties.push_back(MOD::Property("hardest_pT", hardest_jet.pt()));
	   	properties.push_back(MOD::Property("jet_eta", hardest_jet.eta()));
	   	properties.push_back(MOD::Property("pfc_pT", hardest_jet_pfcs[i].pt()));
	   	properties.push_back(MOD::Property("pfc_pdgId", hardest_jet_pfcs[i].user_info<MOD::InfoPFC>().pdgId()));

   		// Now that we've calculated all observables, write them out.

	   string name;
	   
	   int padding = 40;

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

   

   
   
   
}





double angularity_lambda(PseudoJet jet, float k, float beta) {
   
   double lambda = 0.0;

   double R = 0.5;   // Jet Radius.

   double total_pT = 0.0;
   for (unsigned j = 0; j < jet.constituents().size(); j++) {
      total_pT += jet.constituents()[j].pt();
   }

   for (unsigned i = 0; i < jet.constituents().size(); i++) {
      
      PseudoJet constituent = jet.constituents()[i];

      double z_i = constituent.pt() / total_pT;
      
      double delta_R = constituent.delta_R(jet);

      double theta_i = delta_R / R;

      lambda += pow(z_i, k) * pow(theta_i, beta);

   }

   return lambda;

}