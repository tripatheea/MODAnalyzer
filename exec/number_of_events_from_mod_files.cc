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

#include <string>


#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include "../interface/Event.h"

#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>






// using namespace boost::filesystem;

using namespace std;
using namespace fastjet;
using namespace contrib;

using namespace boost::filesystem;


void get_all_files_to_process(std::vector<string> & all_files, boost::filesystem::path input_path);

void get_files_already_processed(std::vector<string> & files_processed, int & total_count, string log_filename);

int main(int argc, char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();


   int number_of_files_to_process = std::numeric_limits<int>::max();
   
   
   if (argc < 3) {
      std::cerr << "You need to give at least two arguments: input_path and output file." << endl;
      return 0;
   }
   else if (argc == 4) {
      cout << "number absent"  << endl;
      number_of_files_to_process = atoi(argv[3]);
   }

   path input_path(argv[1]);
   ofstream output_file(argv[2], ios::out | ios::app);

   cout << endl << endl << "Starting number_of_events_from_mod_files with the following given arguments: " << endl << endl;
   cout << "Input Path   : " << argv[1] << endl;
   cout << "Output File  : " << argv[2] << endl;
   cout << "No. of files : " << number_of_files_to_process << endl;
   
   
   // Recursively collect all filenames to process.

   vector<string> all_filenames;
   get_all_files_to_process(all_filenames, input_path);

   sort(all_filenames.begin(), all_filenames.end());

   // Setup data structures to hold all counts. 
   boost::unordered_map<string, int> event_counts;
   
   int total_files = all_filenames.size();
   int total_validated_events = 0;

   // Read in files we have already counted events for.
   vector<string> files_processed;
   get_files_already_processed(files_processed, total_validated_events, argv[2]);


   int file_counter = 0;
   // Loop through all those files and count events. 
   for (int i = 0; i < total_files; i++) {

      string filename = all_filenames[i];
      filename = filename.substr(filename.length() - 40, 36);

      if ((i >= number_of_files_to_process)) {
         break;
      } 
      else {
         if (std::find(files_processed.begin(), files_processed.end(), filename) == files_processed.end()) {
            
            ifstream file_to_process(all_filenames[i]);
      
            file_counter++;

            if ((file_counter % 1) == 0) {
               cout << "Processing file number " << (file_counter + files_processed.size()) << " / " << total_files << ". I have counted " << total_validated_events << " events so far." << endl;
            }
               
            MOD::Event event_being_read;
            int event_serial_number = 1;

            while (event_being_read.read_event(file_to_process)) {
               
               auto search = event_counts.find(filename);
               if (search != event_counts.end()) {
                  search->second++;
               }
               else {
                  event_counts.insert(make_pair(filename, 1));
               }

               event_being_read = MOD::Event();
               event_serial_number++;
               total_validated_events++;
            }

            output_file << filename << "\t" << event_counts[filename] << endl;   
         }
      }
      
      

      
   }

   
   // vector<string> filenames;
   // for(auto kv : event_counts) {
   //    filenames.push_back(kv.first);
   // }

   // sort( filenames.begin(), filenames.end() );


   // for (string file : filenames) {
   //    output_file << file << "\t" << event_counts[file] << endl;
   // }



   cout << endl << "Everything done." << endl << endl;


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


void get_files_already_processed(std::vector<string> & files_processed, int & total_count, string log_filename) {
   std::ifstream log_file(log_filename);

   std::string line;
   while (std::getline(log_file, line)) {
      std::istringstream iss(line);

      string filename;
      int number;

      iss >> filename >> number;

      files_processed.push_back(filename);

      total_count += number;
   }

}