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

#include "fastjet/ClusterSequence.hh"

#include "../src/event.cc"
#include "../src/fractional_jet_multiplicity.cc"

using namespace std;

bool analyze_event(MODEvent & event_being_read, ofstream & output_file, vector<double> cone_radii, vector<double> pt_cuts);

int main(int argc, char * argv[]) {
   
   int number_of_events_to_process;

   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply three arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory." << std::endl;
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
   ofstream output_file(argv[2], ios::out);

   cout << "Number of outputs to process: " << number_of_events_to_process << endl;
   
   vector<double> cone_radii = {0.3, 0.5, 0.7};
   vector<double> pt_cuts = {50.0, 80.0, 110.0};

   MODEvent event_being_read;

   output_file << "# Event_Number     Run_Number     N_tilde     Jet_Size          Trigger_Name          Fired?     Prescale_1     Prescale_2     Cone_Radius     pT_Cut     Hardest_pT" << endl;

   int event_serial_number = 1;
   while( ( event_being_read.read_event(data_file, event_being_read) ) && ( event_serial_number <= number_of_events_to_process ) ) {
      
      cout << "Processing event number " << event_serial_number << endl;

      analyze_event(event_being_read, output_file, cone_radii, pt_cuts);
      event_serial_number++;
   }

   return 0;
}


bool analyze_event(MODEvent & event_being_read, ofstream & output_file, vector<double> cone_radii, vector<double> pt_cuts) {

   // Retrieve the assigned trigger and store information about that trigger (prescales, fired or not).
   // Also calculate everything and record those along with the trigger information.

   string assigned_trigger_name = event_being_read.assigned_trigger_name();
   const MODTrigger assigned_trigger = event_being_read.trigger_by_name(assigned_trigger_name);

   pair<int, int> prescales = assigned_trigger.prescale_pair();
   
   bool fired = assigned_trigger.fired();
   

   int prescale_1 = prescales.first;
   int prescale_2 = prescales.second;

   double hardest_pt = event_being_read.hardest_pt();

   // Calculate everything for each value of R and pt_cut.

   for(unsigned int r = 0; r < cone_radii.size(); r++) {
      for (unsigned int p = 0; p < pt_cuts.size(); p++) {

         // Calculate N_tilde.
         FractionalJetMultiplicity ntilde = FractionalJetMultiplicity(cone_radii[r], pt_cuts[p]);
         double N_tilde = ntilde(event_being_read.pseudojets());
         
         // Calculate jet size (fastjet)
         JetDefinition jet_def(antikt_algorithm, cone_radii[r]);
         ClusterSequence cs(event_being_read.pseudojets(), jet_def);
         vector<PseudoJet> jets = cs.inclusive_jets(pt_cuts[p]);

         output_file << setw(12) << event_being_read.event_number()
                  << setw(15) << event_being_read.run_number()
                  << setw(14) << showpoint << setprecision(6) << N_tilde
                  << setw(9) << jets.size()
                  << setw(35) << assigned_trigger_name
                  << setw(5) << fired
                  << setw(12) << prescale_1
                  << setw(15) << prescale_2
                  << setw(18) << setprecision(2) << cone_radii[r]
                  << setw(12) << noshowpoint << setprecision(3) << pt_cuts[p]
                  << setw(16) << showpoint << setprecision(8) << hardest_pt
                  << endl;             
      }
   }
}