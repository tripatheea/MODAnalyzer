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
      
      if( (event_serial_number % 5000) == 0 )
         cout << "Processing event number " << event_serial_number << endl;

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

   fastjet::PseudoJet trigger_jet = event_being_read.trigger_jet();
   fastjet::PseudoJet closest_jet_to_trigger_jet = event_being_read.closest_fastjet_jet_to_trigger_jet();

   vector<MOD::Property> properties;
   
   if ( (event_being_read.cms_jets().size() == 0) or (event_being_read.jets().size() == 0) or ( ! trigger_jet.has_user_info())) {
      return;
   }

   vector<MOD::Trigger> triggers = event_being_read.triggers();
  

   
   try {
      for (unsigned i = 0; i < triggers.size(); i++) {
         
         // if (triggers[i].fired()) {

            fastjet::PseudoJet trigger_jet = event_being_read.trigger_jet();

            properties.push_back(MOD::Property("# Entry", "  Entry"));

            properties.push_back(MOD::Property("Event_Number", event_being_read.event_number()));
            properties.push_back(MOD::Property("Run_Number", event_being_read.run_number()));

            properties.push_back(MOD::Property("trig_jet_matched", (int) event_being_read.trigger_jet_is_matched())); 
            properties.push_back(MOD::Property("jet_quality", trigger_jet.user_info<MOD::InfoCalibratedJet>().jet_quality())); 
   
            properties.push_back(MOD::Property("Cor_Hardest_pT", trigger_jet.pt() * trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC()));  
            properties.push_back(MOD::Property("hardest_eta", trigger_jet.eta()));  
            properties.push_back(MOD::Property("Prescale", triggers[i].prescale()));
            properties.push_back(MOD::Property("Trigger_Name", triggers[i].name()));         
       
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

            properties.clear();

         // }
      }
   }
   catch (exception& e) {

   }   


   
   

   


   

   
}