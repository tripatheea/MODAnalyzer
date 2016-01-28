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

void analyze_luminosity(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, set<int> & lumi_blocks);

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

   set<int> lumi_blocks;

   MOD::Event event_being_read;

   int event_serial_number = 1;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {
      
      if( (event_serial_number % 5000) == 0 )
         cout << "Processing event number " << event_serial_number << endl;

      analyze_luminosity(event_being_read, output_file, event_serial_number, lumi_blocks);
      
      event_being_read = MOD::Event();
      event_serial_number++;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_luminosity(MOD::Event & event_being_read, ofstream & output_file, int & event_serial_number, set<int> & lumi_blocks) {


   int lumi_block = event_being_read.condition().lumi_block();

   const bool lumi_block_already_appeared = lumi_blocks.find(lumi_block) != lumi_blocks.end();

   if ( ! lumi_block_already_appeared) {

      lumi_blocks.insert(lumi_block);

      vector<MOD::Property> properties;

      properties.push_back(MOD::Property("# Entry", "  Entry"));

      properties.push_back(MOD::Property("Event_Number", event_being_read.event_number()));
      properties.push_back(MOD::Property("Run_Number", event_being_read.run_number()));

      properties.push_back(MOD::Property("lumi_block", lumi_block)); 
      properties.push_back(MOD::Property("valid_lumi", event_being_read.condition().valid_lumi())); 

      properties.push_back(MOD::Property("avg_inst_lumi", event_being_read.condition().average_instantaneous_lumi())); 
      properties.push_back(MOD::Property("intg_del_lumi", event_being_read.condition().integrated_delivered_lumi())); 
      properties.push_back(MOD::Property("intg_rec_lumi", event_being_read.condition().integrated_recorded_lumi())); 
      properties.push_back(MOD::Property("time", event_being_read.condition().time())); 

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

      properties.clear();
   }


}