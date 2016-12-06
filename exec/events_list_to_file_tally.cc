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

string find_correct_file_for_event(int run_number, int event_number, unordered_map<string, string> & registry_info);
void load_registry_information(string registry_filename, unordered_map<string, string> & registry_info);

int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();


   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply five arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }


   string input_log(argv[1]);
   string output_log(argv[2]);
   string registry_file(argv[3]);

   
   cout << endl << endl << "Starting process with the following given arguments: " << endl;
   cout << "Input Log    : " << argv[1] << endl;
   cout << "Output Log   : " << argv[2] << endl;
   cout << "Registry File: " << argv[3] << endl;

   unordered_map<string, string> registry_info;


   load_registry_information(registry_file, registry_info);

   cout << "Finished loading registry to memory." << endl;

   unordered_map<string, int> files_events_count;

   ifstream input_file(input_log);
   ofstream output_file(output_log, ios::out);

   int run_number, event_number;
   string line;
   while(getline(input_file, line)) {
      stringstream iss(line);

      iss >> run_number >> event_number;

      string filename = find_correct_file_for_event(run_number, event_number, registry_info);
      auto search = files_events_count.find(filename);
      
      if (search != files_events_count.end()) {
         search->second++;
      }
      else {
         files_events_count.insert(make_pair(filename, 1));
      }

   }


   for (auto kv : files_events_count) {
      output_file << kv.first << "\t" << kv.second << endl;
   }

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}



string find_correct_file_for_event(int run_number, int event_number, unordered_map<string, string> & registry_info) {

   auto search = registry_info.find(to_string(run_number) + "_" + to_string(event_number));
   
   if(search != registry_info.end()) {
      // std::cout << "Found " << search->first << " " << search->second << '\n';
      return search->second;
   }
   else {
      // cout << "Oops didn't find it." << endl;
   }

   return "";
}





void load_registry_information(string registry_filename, unordered_map<string, string> & registry_info) {

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

      registry_info.emplace(to_string(registry_run_number) + "_" + to_string(registry_event_number), root_filename.substr(0, 36));  // 36 because the filenames without extensions are 37 characters long.

      line_number++;
   }
}