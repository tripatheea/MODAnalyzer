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

#include "../interface/Event.h"
#include "../interface/Property.h"

using namespace std;

void analyze_event(MOD::Event & event_being_read, string output_path, unordered_map<int, string> & registry_info, ofstream & completed_events_file_output, boost::unordered_map<pair<int, int>, int> & completed_events);
string find_correct_file_for_event(int event_number, unordered_map<int, string> & registry_info);
void load_registry_information(string registry_filename, unordered_map<int, string> & registry_info);


int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();


   if (argc <= 3) {
        std::cerr << "ERROR: You need to supply five arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }


   string input_path(argv[1]);
   string current_root_file(argv[2]);
   string registry_file(argv[3]);
   string output_path(argv[4]);
   string completed_log_filename(argv[5]);

   ifstream input_mod_file(input_path + "/" + current_root_file);
   
   
   cout << endl << endl << "Starting process with the following given arguments: " << endl;
   cout << "Input Path: " << argv[1] << endl;
   cout << "Input Filename: " << input_path << "/" << argv[2] << endl;
   cout << "Registry File : " << argv[3] << endl;
   cout << "Output Path   : " << argv[4] << endl;
   cout << "Completed Log : " << argv[5] << endl << endl << endl;

   unordered_map<int, string> registry_info;


   load_registry_information(registry_file, registry_info);

   cout << "Finished loading registry to memory." << endl;


   boost::unordered_map<pair<int, int>, int> completed_events;

   // Load completed events to the vector "completed_events_"
   ifstream completed_file(completed_log_filename.c_str());
   
   string line;
   int line_number = 1;
   while(getline(completed_file, line)) {
      
      if (line_number % 100000 == 0)
         cout << "On line number " << line_number << endl;

      istringstream iss(line);
      int event_number, run_number;
      iss >> run_number >> event_number;
      
      completed_events.insert(make_pair(make_pair(run_number, event_number), 1));

      line_number++;
   }

   cout << "Finished loading \"completed\" events to memory." << endl;

   ofstream completed_events_file_output;
   completed_events_file_output.open(completed_log_filename.c_str(), ios::out | ios::app );


   MOD::Event event_being_read;

   int event_serial_number = 1;
   while ( event_being_read.read_event(input_mod_file) ) {
      
      if( (event_serial_number % 5000) == 0 )
         cout << "Processing event number " << event_serial_number << endl;

      
      analyze_event(event_being_read, output_path, registry_info, completed_events_file_output, completed_events);
      
      event_being_read = MOD::Event();
      event_serial_number++;

   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}


void analyze_event(MOD::Event & event_being_read, string output_path, unordered_map<int, string> & registry_info, ofstream & completed_events_file_output, boost::unordered_map<pair<int, int>, int> & completed_events) {

   ofstream log_file(output_path + "/error.log", ios::out | ios::app);

   if (completed_events.find(make_pair(event_being_read.run_number(), event_being_read.event_number())) == completed_events.end()) {
      
      // Event not processed already.

      string correct_filename = find_correct_file_for_event(event_being_read.event_number(), registry_info);

      if (correct_filename == "") {
         log_file << "Could not find the correct filename for the following event: " << event_being_read.event_number() << " " << event_being_read.run_number() << endl;
         return;
      }
      else {
         ofstream output_file(output_path + "/" + correct_filename + ".mod", ios::out | ios::app);

         // See if the event has already been added to the list. 
         
         // cout << event_being_read << endl;
         // cout << "Writing file " << output_path << "/" << correct_filename << ".mod" << endl;
         output_file << event_being_read;
   
      }

      // Add it to the map completed_events and to the file.
      completed_events.insert(make_pair(make_pair(event_being_read.run_number(), event_being_read.event_number()), 1));

      // completed_events_file_output << event_being_read.run_number() << "\t" << event_being_read.event_number() << endl; 

   }
   else {
      // cout << "Event " << event_being_read.run_number() << "\t" << event_being_read.event_number() << " already processed so skipping." << endl;
   }
   
}

string find_correct_file_for_event(int event_number, unordered_map<int, string> & registry_info) {

   auto search = registry_info.find(event_number);
   
   if(search != registry_info.end()) {
      // std::cout << "Found " << search->first << " " << search->second << '\n';
      return search->second;
   }
   else {
      //cout << "Oops didn't find it." << endl;
   }

   return "";
}

void load_registry_information(string registry_filename, unordered_map<int, string> & registry_info) {

   cout << "Loading registry to memory." << endl;

   ifstream registry(registry_filename);

   int line_number = 1;

   string line;
   // while ((getline(registry, line)) && (line_number < 100000000)) {
   while (getline(registry, line)) {

      if (line_number % 100000 == 0)
         cout << "On line number " << line_number << endl;

      istringstream iss(line);

      int registry_event_number, registry_run_number;
      string file_path, root_filename;

      iss >> registry_event_number >> registry_run_number >> file_path >> root_filename;

      registry_info.emplace(registry_event_number, root_filename.substr(0, 36));  // 36 because the filenames without extensions are 37 characters long.

      line_number++;
   }
}