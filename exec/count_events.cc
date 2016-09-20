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
#include "../interface/Event.h"

#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>


// using namespace boost::filesystem;

using namespace std;
using namespace fastjet;

using namespace boost::filesystem;


void count_events(MOD::Event & event_being_read, boost::unordered_map<int, int> & npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut);
void get_all_files_to_process(std::vector<string> & all_files, boost::filesystem::path input_path);

int main(int argc, char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();

   int number_of_events_to_process = std::numeric_limits<int>::max();
   int number_of_files_to_process = std::numeric_limits<int>::max();

   path input_path(argv[1]);
   bool use_latex(argv[2]);

   
   if (argc == 4) {
      // Number of files to process given. 
      number_of_files_to_process = atoi(argv[3]);
   }
   else if (argc == 5) {
      number_of_files_to_process = atoi(argv[3]);
      number_of_events_to_process = atoi(argv[4]);
   }
   else if (argc < 3) {
      std::cerr << "You need to give at least two arguments: input_path and whether or not to use LaTeX." << endl;
      return 0;
   }

   cout << use_latex << endl;
   
   cout << endl << endl << "Starting count_events with the following given arguments: " << endl << endl;
   cout << "Input Path   : " << argv[1] << endl;
   cout << "Use LaTeX    : " << use_latex << endl;
   cout << "No. of Files : " << number_of_files_to_process << endl;
   cout << "No. of Events: " << number_of_events_to_process << endl;
   
   // Recursively collect all filenames to process.

   vector<string> all_filenames;
   get_all_files_to_process(all_filenames, input_path);

   // Setup data structures to hold all counts. 
   boost::unordered_map<int, int> npv;
   boost::unordered_map<string, int> jet_quality;
   boost::unordered_map<float, int> passes_softdrop;

   // int total_files = all_filenames.size();
   int total_files = all_filenames.size();
   int total_validated_events = 0;
   int assigned_trigger_fired = 0;
   int ak5_match = 0;
   int eta_24_cut = 0;
   

   int file_counter = 0;
   // Loop through all those files and count events. 
   for (int i = 0; i < total_files; i++) {

      if (i < number_of_files_to_process) {
         ifstream file_to_process(all_filenames[i]);

         file_counter++;

         if ((file_counter % 100) == 0)
            cout << "Processing file number " << file_counter << " / " << total_files << endl;

         MOD::Event event_being_read;
         int event_serial_number = 1;

         while ((event_being_read.read_event(file_to_process)) && (event_serial_number <= number_of_events_to_process)) {
         // while (event_being_read.read_event(file_to_process)) {

            // if ((event_serial_number % 1000) == 0)
               cout << "Reading event number " << event_serial_number << endl;
            
            count_events(event_being_read, npv, jet_quality, passes_softdrop, assigned_trigger_fired, ak5_match, eta_24_cut);

            event_being_read = MOD::Event();
            event_serial_number++;
            total_validated_events++;
         }   
      }
      else {
         break;
      }
      
   }

   
   cout << endl << "Everything done. Printing summary." << endl << endl;

   cout << endl << "==================================================================" << endl << endl;

   int grand_total = 20022826;


   cout << "Assigned Trigger Fired: " << assigned_trigger_fired << endl << endl;
   cout << "AK5 Match: " << ak5_match << endl << endl;
   cout << "|eta| < 2.4: " << eta_24_cut << endl << endl;

   cout << "====== NPV ======" << endl;
   for (auto kv : npv) {
      cout << kv.first << " => " << kv.second << endl;
   }
   cout << endl << endl;

   cout << "====== Jet Quality ======" << endl;
   for (auto kv : jet_quality) {
      cout << kv.first << " => " << kv.second << endl;
   }
   cout << endl << endl;

   cout << "====== Total ======" << endl;
   cout << "Total Validated Events: " << total_validated_events << endl;   


   if (use_latex) {
      cout << "\\hline" << endl;
      cout << "\\hline" << endl;
      cout << "& Events & Fraction\\\\" << endl;
      cout << "\\hline" << endl;
      cout << "Jet Primary Dataset & 20,022,826 & 1.0 \\\\" << endl;
      cout << "Validated Run & " << total_validated_events << " & " << fixed << setw(5) << (total_validated_events/grand_total) << " \\\\" << endl;
      cout << "Loose Jet Quality (\\Tab{tab:jet_quality}) & " << jet_quality["LOOSE"] << " & \at{?} \\\\" << endl;
      cout << "Assigned Trigger Fired (\\Tab{tab:trigger_table}) & 444,823 & 0.022 \\\\" << endl;
      cout << "\\hline" << endl;
      cout << "AK5 Match & \at{?} & \at{?}\\\\" << endl;
      cout << "$\\left| \eta \right| < 2.4$ & \at{?} & \at{?} \\\\" << endl;
      cout << "\\hline" << endl;
      cout << "Passes Soft Drop ($z_g > z_{\rm cut}$) & \at{?} & \at{?} \\\\" << endl;
      cout << "\\hline" << endl;
      cout << "\\hline" << endl;
   }
   

}


void get_all_files_to_process(std::vector<string> & all_files, boost::filesystem::path input_path) {
   
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
         get_all_files_to_process(all_files, itr->path());
      }
   }

}

void count_events(MOD::Event & event_being_read, boost::unordered_map<int, int> & npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut) {
   
      
   // NPV
   auto search = npv.find(event_being_read.condition().npv());

   if (search != npv.end())
      search->second++;
   else
      npv.insert(make_pair(event_being_read.condition().npv(), 1));  

   // First, find the trigger jet.
   PseudoJet trigger_jet = event_being_read.trigger_jet();

   // Jet Quality
   if (trigger_jet.has_user_info()) {
      int quality = trigger_jet.user_info<MOD::InfoCalibratedJet>().jet_quality();
      string quality_string;

      if (quality == 0)
         quality_string = "FAILED";
      else if (quality == 1)
         quality_string = "LOOSE";
      else if (quality == 2)
         quality_string = "MEDIUM";
      else if (quality == 3)
         quality_string = "TIGHT";
      
      auto search = jet_quality.find(quality_string);

      if (search != jet_quality.end()) {
         
         search->second++;

         if (quality_string == "TIGHT") {
            jet_quality.find("MEDIUM")->second++;
            jet_quality.find("LOOSE")->second++;
         }
         else if (quality_string == "MEDIUM") {
            jet_quality
         }

      }
      else {
         jet_quality.insert(make_pair(quality_string, 1));
      }
   }
   


   // Passes SoftDrop

   if (event_being_read.assigned_trigger_fired())
      assigned_trigger_fired++;

   if (event_being_read.trigger_jet_is_matched())
      ak5_match++;

   if (abs(event_being_read.hardest_jet().eta()) < 2.4)
      eta_24_cut++;

   // PseudoJet trigger_jet = event_being_read.trigger_jet();
   // event_being_read.convert_to_one_jet();


   // if (event_being_read.trigger_jet_is_matched() && (trigger_jet.user_info<MOD::InfoCalibratedJet>().jet_quality() >= 1)) {   // Jet quality level: FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3      
   //    output_file << event_being_read;
   // }
   
}
