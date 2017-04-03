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

#include <boost/unordered_map.hpp>


#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>


#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"


#include "../interface/Event.h"
#include "../interface/Property.h"

using namespace std;
using namespace fastjet;
using namespace contrib;
using namespace boost::filesystem;


void analyze_pfc(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number);

void get_all_files_to_process(std::vector<string> & all_files, vector<boost::filesystem::path> input_paths);

double angularity_lambda(PseudoJet jet, float k, float beta);

int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();

   ofstream output_file(argv[1], ios::out | ios::app);

   
 
   if (argc < 3) {
      std::cerr << "You need to give at least two arguments: output path and one or more input directories." << endl;
      return 0;
   }

   vector<path> input_paths;
   for (int i=2; i < argc; i++) {
     input_paths.push_back(argv[i]);
   }


   
   cout << endl << endl << "Starting count_events with the following given arguments: " << endl << endl;
   cout << "Output File     : " << argv[1] << endl;
   cout << "Input Paths (" << input_paths.size() << ") : ";
   for (unsigned i=0; i < input_paths.size(); i++) {
    if (i > 0)
        cout << "                  " << input_paths[i] << endl;
    else 
        cout << input_paths[i] << endl;
   }   

   // Recursively collect all filenames to process.

   vector<string> all_filenames;
   get_all_files_to_process(all_filenames, input_paths);

   // Sort the list.
   std::sort(all_filenames.begin(), all_filenames.end());
   
   cout << endl << endl << "Starting counting event on " << all_filenames.size() << " files." << endl << endl;
   

   int file_counter = 0;
   int event_counter = 0;
   // Loop through all those files and count events. 
   for (unsigned i = 0; i < all_filenames.size(); i++) {

		ifstream file_to_process(all_filenames[i]);

		file_counter++;

		if ((file_counter % 100) == 0)
			cout << "Processing file number " << file_counter << " / " << all_filenames.size() << endl;


	   MOD::Event event_being_read;

	   int event_serial_number = 1;
	   while( event_being_read.read_event(file_to_process) ) {
	      
	      if( (event_counter % 10000) == 0 )
	         cout << "Processed " << (event_counter - 1) << " events so far." << endl;

	      // Write out version info in the output file for the "syncing plots" thing to work (as it needs to figure out which directory to put things into).
	      if (event_serial_number == 1)
	         output_file << "%" << " Version " << event_being_read.version() << endl;

	    
	      analyze_pfc(event_being_read, output_file, event_serial_number);
	      
	      event_being_read = MOD::Event();
	      event_serial_number++;
	      event_counter++;

	   }

	   event_being_read = MOD::Event();
	   event_serial_number++;
	// }
  
   }


   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_counter - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}









void get_all_files_to_process(std::vector<string> & all_files, vector<boost::filesystem::path> input_paths) {
   
   for (unsigned i = 0; i < input_paths.size(); i++) {
       boost::filesystem::path input_path = input_paths[i];
       
       directory_iterator end_itr;

       for (directory_iterator itr(input_path); itr != end_itr; ++itr) {
          
          if (is_regular_file(itr->path())) {
             string current_file = itr->path().string();
             
             if (current_file.substr( current_file.length() - 3, current_file.length()) == "mod") {
                all_files.push_back(current_file);   
             }

          }
          else {
             // cout << itr->path().string() << endl;
             get_all_files_to_process(all_files, { itr->path() });
          }
       }
    }
}

void analyze_pfc(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number) {


   JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
   PseudoJet hardest_jet = event_being_read.hardest_jet();

   if ( ! (hardest_jet.E() > 0.0)) {
      return;
   }

   double jec = event_being_read.get_hardest_jet_jec();
   
 
   vector<PseudoJet> hardest_jet_pfcs = hardest_jet.constituents();
      
   for (unsigned i = 0; i < hardest_jet_pfcs.size(); i++) {

   		vector<MOD::Property> properties;

   		properties.push_back(MOD::Property("# Entry", "  Entry"));

        properties.push_back(MOD::Property("event_number", event_being_read.event_number()));
        properties.push_back(MOD::Property("prescale", event_being_read.weight()));
	   	properties.push_back(MOD::Property("hardest_pT", jec * hardest_jet.pt()));
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