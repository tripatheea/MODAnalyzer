#include <iostream>
#include <unordered_map>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>
#include <iomanip>
#include <chrono>

#include <cstdio>

#include "fastjet/ClusterSequence.hh"
#include "../interface/event.h"
#include "../interface/fractional_jet_multiplicity.h"

using namespace std;
using namespace fastjet;

void skim(MOD::Event & event_being_read, ofstream & output_file);
bool jets_match(MOD::Event & event_being_read);
bool pseudojets_compare(PseudoJet a, PseudoJet b);

int main(int argc, char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();

   int number_of_events_to_process;

   if (argc <= 3) {
        std::cerr << "ERROR: You need to supply four arguments- first, path to the input data; second, path to the output file; ; third, path to save the log file to, fourth, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./filter (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }
   else if (argc == 4) {
      // Fourth argument is missing, process everything.
      number_of_events_to_process = std::numeric_limits<int>::max();
   }
   else {
      // Fourth argument gives the number of events to process.
      number_of_events_to_process = stoi(argv[4]);
   }

   ofstream log_stream(argv[3], ios::out | ios::app );
   ofstream num_stream( std::string(argv[3]) + ".num", ios::out | ios::app );

   // connect stream buffers
   std::streambuf *cerrbuf = std::cerr.rdbuf();
   std::cerr.rdbuf(log_stream.rdbuf () );

   
   ifstream data_file(argv[1]);
   ofstream output_file(argv[2], ios::out | ios::app );
   
   cout << endl << endl << "Starting skimming with the following given arguments: " << endl;
   cout << "Input file: " << argv[1] << endl;
   cout << "Output file: " << argv[2] << endl;
   cout << "Number of events: ";
   if(argc == 4)
      cout << "ALL" << endl << endl;
   else
      cout << number_of_events_to_process << endl << endl;

   MOD::Event event_being_read;

   int event_serial_number = 1;
   int events_with_mismatched_jets = 0;
   while( event_being_read.read_event(data_file) && ( event_serial_number <= number_of_events_to_process ) ) {

      if( (event_serial_number % 1000) == 0 ) {
         cout << "Skimming event number " << event_serial_number << endl;
      }

      skim(event_being_read, output_file);
      
      if( ! jets_match(event_being_read)) {
         events_with_mismatched_jets++;
         cerr << "AK5 Jets don't match for event number: " << event_being_read.event_number() << " run number: " << event_being_read.run_number() << " in file: " << argv[1] << endl;
      }

      event_being_read = MOD::Event();
      event_serial_number++;
   }

   

   // cout << events_with_mismatched_jets << " events have mismatched AK5 jets!" << endl;

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();

   cout << "Finished skimming " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl << endl;

   // Write number of events processed to the num stream so that we can sum it up later.
   num_stream << (event_serial_number - 1) << endl;

   // restore
   std::cout.flush ();
   std::cout.rdbuf(cerrbuf);

}


void skim(MOD::Event & event_being_read, ofstream & output_file) {
   
   if (event_being_read.assigned_trigger_fired())
      output_file << event_being_read.make_string();
}

bool jets_match(MOD::Event & event_being_read) {
   
   double pt_cut = 3.00;
   double cone_radius = 0.5;
   
   vector<PseudoJet> cms_jets = event_being_read.CMS_pseudojets();
   vector<PseudoJet> pfcandidates = event_being_read.pseudojets();

   // Cluster the pfcandidates using Fastjet.
   JetDefinition jet_def(antikt_algorithm, cone_radius);
   ClusterSequence cs(pfcandidates, jet_def);
   vector<PseudoJet> fastjet_jets = cs.inclusive_jets(pt_cut);

   // Compare the number of jets first.
   if (cms_jets.size() != fastjet_jets.size()) {
      cerr << endl << "Different jet size for CMS vs. FastJet; " << cms_jets.size() << " vs " << fastjet_jets.size() << "." << endl;
      return false;
   }


   // Next, compare if fastjet_jets with cms_jets or not upto 10e-4 precision.
   
   unsigned max_number_of_jets =  max(fastjet_jets.size(), cms_jets.size());
  
   // First sort both vectors by px so that we can compare them one by one.
   sort(fastjet_jets.begin(), fastjet_jets.end(), pseudojets_compare);
   sort(cms_jets.begin(), cms_jets.end(), pseudojets_compare);

   double tolerance = pow(10, -3);
   for (unsigned i = 0; i < max_number_of_jets; i++) {
      
      // Next check if CMS AK5 jets match locally clustered FastJet AK5 jets within a certain precision. 
      
      if ( ( abs(fastjet_jets[i].px() - cms_jets[i].px()) > tolerance ) || ( abs(fastjet_jets[i].py() - cms_jets[i].py()) > tolerance ) || ( abs(fastjet_jets[i].pz() - cms_jets[i].pz()) > tolerance ) || ( abs(fastjet_jets[i].E() - cms_jets[i].E()) > tolerance ) ) {
         
         cerr << endl << fixed << setprecision(5) << "FastJet: " <<  fastjet_jets[i].px() << "   " << fastjet_jets[i].py() << "   " << fastjet_jets[i].pz() << "   " << fastjet_jets[i].E() << endl;
         cerr         << fixed << setprecision(5) << "CMS:     " << cms_jets[i].px() << "   " << cms_jets[i].py() << "   " << cms_jets[i].pz() << "   " << cms_jets[i].E() << endl;
         
         // Jet quality level set to blank and eta_cut set to 0.0 because they're not going to be used anyway plus this will break the output if there's an error in the code somewhere making it easier to debug.
         MOD::CalibratedJet hardest_uncorrected_jet = event_being_read.hardest_jet(false, false, false, "", 0.0);

         if( ( abs(hardest_uncorrected_jet.pseudojet().px() - cms_jets[i].px()) < tolerance ) || ( abs(hardest_uncorrected_jet.pseudojet().py() - cms_jets[i].py()) < tolerance ) || ( abs(hardest_uncorrected_jet.pseudojet().pz() - cms_jets[i].pz()) < tolerance ) || ( abs(hardest_uncorrected_jet.pseudojet().E() - cms_jets[i].E()) < tolerance ) ) {
            cerr << "This error is with the hardest jet!" << endl;
         }

         return false;

      }
   }

   return true;  
}

bool pseudojets_compare(PseudoJet a, PseudoJet b) {
   if (a.pt() > b.pt())
      return true;
   return false;
}