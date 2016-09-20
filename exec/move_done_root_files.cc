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


   string root_files_path(argv[1]);
   string done_root_files_path(argv[2]);
   string registry_file(argv[3]);
   string completed_log_filename(argv[4]);

   // ifstream input_mod_file(input_path + "/" + current_root_file);
   
   
   cout << endl << endl << "Starting process with the following given arguments: " << endl;
   cout << "Root Files Path : " << argv[1] << endl;
   cout << "Done Root Files : " << argv[2] << endl;
   cout << "Registry File   : " << argv[3] << endl;
   cout << "Completed Log : " << argv[4] << endl << endl << endl;

   unordered_map<int, string> registry_info;

   load_registry_information(registry_file, registry_info);

   cout << "Finished loading registry to memory." << endl;

   boost::unordered_map<string, int> done_root_files;


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
      
      // done_root_files.push_back(find_correct_file_for_event(event_number, registry_info));
      done_root_files.insert(make_pair(find_correct_file_for_event(event_number, registry_info), 1));

      line_number++;
   }

   cout << "Finished loading \"completed\" events to memory." << endl;
   cout << "Moving them now." << endl;

   vector<string> unique_done_filenames; 

   for (auto kv : done_root_files) {
      unique_done_filenames.push_back(kv.first);
   } 

   cout << unique_done_filenames.size() << endl;
   
   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing in " << elapsed_seconds << " seconds!" << endl;

   return 0;
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