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
#include <utility>
#include <chrono>
#include <algorithm>

#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>

#include "../interface/Event.h"
#include "../interface/Property.h"

using namespace std;
using namespace boost::filesystem;

void analyze_event(MOD::Event & event_being_read, ofstream & completed_events_file_output, boost::unordered_map<string, int> & completed_events, int & number_of_events_written);

void get_all_files_to_process(std::vector<string> & all_files, boost::filesystem::path input_path);

int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();


   if (argc <= 1) {
        std::cerr << "ERROR: You need to supply five arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }


   string input_path(argv[1]);
   string output_file(argv[2]);
   

   
   cout << endl << endl << "Starting process with the following given arguments: " << endl;
   cout << "Input Path: " << argv[1] << endl;
   cout << "Output File   : " << argv[2] << endl;

   boost::unordered_map<string, int> completed_events;

   // Load completed events to the vector "completed_events_"
   ifstream completed_file(output_file.c_str());
   
   string line;
   int line_number = 1;
   while(getline(completed_file, line)) {
      
      if (line_number % 100000 == 0)
         cout << "On line number " << line_number << endl;

      istringstream iss(line);
      int event_number, run_number;
      iss >> run_number >> event_number;
      
      completed_events.insert(make_pair(to_string(run_number) + "_" + to_string(event_number), 1));

      line_number++;
   }

   cout << "Finished loading \"completed\" events to memory." << endl;

   ofstream completed_events_file_output;
   completed_events_file_output.open(output_file.c_str(), ios::out | ios::app );


   vector<string> all_filenames;
   get_all_files_to_process(all_filenames, input_path);

   int file_counter = 0;
   for (unsigned i=0; i < all_filenames.size(); i++) {
      ifstream input_mod_file(all_filenames[i]);

      file_counter++;

      // if ((file_counter % 100) == 0)
      //    cout << "Processing file number " << file_counter << " / " << all_filenames.size() << endl;

      cout << "Processing file " << (i + 1) << " / " << all_filenames.size() << ": " << all_filenames[i] << endl;
      
      MOD::Event event_being_read;

      int event_serial_number = 1;
      int number_of_events_written = 0;
      while ( event_being_read.read_event(input_mod_file) ) {
         
         // if( (event_serial_number % 5000) == 0 )
         //    cout << "Processing event number " << event_serial_number << endl;

         
         analyze_event(event_being_read, completed_events_file_output, completed_events, number_of_events_written);
         
         event_being_read = MOD::Event();
         event_serial_number++;
      }

   }


   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << all_filenames.size() << " files in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_event(MOD::Event & event_being_read, ofstream & completed_events_file_output, boost::unordered_map<string, int> & completed_events, int & number_of_events_written) {


   if (completed_events.find(to_string(event_being_read.run_number()) + "_" + to_string(event_being_read.event_number())) == completed_events.end()) {
      
      // Event not processed already.

      completed_events_file_output << event_being_read.run_number() << "\t" << event_being_read.event_number() << endl;

      number_of_events_written++;

      // Add it to the map completed_events and to the file.
      completed_events.insert(make_pair(to_string(event_being_read.run_number()) + "_" + to_string(event_being_read.event_number()), 1));

   }
   else {
      cout << "Event " << event_being_read.run_number() << "\t" << event_being_read.event_number() << " already processed so skipping." << endl;
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

