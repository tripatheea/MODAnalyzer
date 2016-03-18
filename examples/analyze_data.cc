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

void analyze_event_data(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number);

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
      
      if( (event_serial_number % 100) == 0 )
         cout << "Processing event number " << event_serial_number << endl;

      // Write out version info in the output file for the "syncing plots" thing to work (as it needs to figure out which directory to put things into).
      if (event_serial_number == 1)
         output_file << "%" << " Version " << event_being_read.version() << endl;

      if (event_being_read.assigned_trigger_fired())
         analyze_event_data(event_being_read, output_file, event_serial_number);
      
      event_being_read = MOD::Event();
      event_serial_number++;

   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_event_data(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number) {

	fastjet::PseudoJet trigger_jet = event_being_read.trigger_jet();
	fastjet::PseudoJet closest_jet_to_trigger_jet = event_being_read.closest_fastjet_jet_to_trigger_jet();

	vector<MOD::Property> properties;

	properties.push_back(MOD::Property("# Entry", "  Entry"));

	properties.push_back(MOD::Property("event_number", event_being_read.event_number()));
	properties.push_back(MOD::Property("run_number", event_being_read.run_number()));

	properties.push_back(MOD::Property("trig_jet_matched", (int) event_being_read.trigger_jet_is_matched())); 
	properties.push_back(MOD::Property("jet_quality", trigger_jet.user_info<MOD::InfoCalibratedJet>().jet_quality())); 

	properties.push_back(MOD::Property("uncor_hardest_pT", trigger_jet.pt()));
	properties.push_back(MOD::Property("cor_hardest_pT", trigger_jet.pt() * trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC()));

	properties.push_back(MOD::Property("prescale", event_being_read.assigned_trigger_prescale()));
	properties.push_back(MOD::Property("trigger_name", event_being_read.assigned_trigger_name()));

	properties.push_back( MOD::Property("no_of_const", trigger_jet.user_info<MOD::InfoCalibratedJet>().number_of_constituents()) );
	properties.push_back( MOD::Property("chrg_multip", trigger_jet.user_info<MOD::InfoCalibratedJet>().charged_multiplicity()) );
	properties.push_back( MOD::Property("neu_had_frac", trigger_jet.user_info<MOD::InfoCalibratedJet>().neutral_hadron_fraction()) );
	properties.push_back( MOD::Property("neu_em_frac", trigger_jet.user_info<MOD::InfoCalibratedJet>().neutral_em_fraction()) );
	properties.push_back( MOD::Property("chrg_had_frac", trigger_jet.user_info<MOD::InfoCalibratedJet>().charged_hadron_fraction()) );
	properties.push_back( MOD::Property("chrg_em_frac", trigger_jet.user_info<MOD::InfoCalibratedJet>().charged_em_fraction()) );
	
	properties.push_back( MOD::Property("JEC", trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC()) );

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