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


#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>


#include "../interface/Event.h"
#include "../interface/Property.h"

using namespace std;
using namespace boost::filesystem;


void get_all_files_to_process(std::vector<string> & all_files, boost::filesystem::path input_path);
void load_registry_information(string registry_filename, boost::unordered_map<string, string> & registry_info);

string find_correct_file_for_event(int run_number, int event_number, boost::unordered_map<string, string> & registry_info);


int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();


   if (argc <= 3) {
        std::cerr << "ERROR: You need to supply five arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }


   string input_path(argv[1]);
   string output_path(argv[2]);
   string registry_file(argv[3]);
   string completed_log_filename(argv[4]);
 
   
   cout << endl << endl << "Starting process with the following given arguments: " << endl;
   cout << "Input Path    : " << argv[1] << endl;
   cout << "Output Path   : " << argv[2] << endl;
   cout << "Registry File : " << argv[3] << endl;
   cout << "Completed Log : " << argv[4] << endl << endl << endl;


   boost::unordered_map<string, string> registry_info;
   load_registry_information(registry_file, registry_info);


   // boost::unordered_map<string, int> completed_events;

   // Load completed events to the vector "completed_events_"
   // ifstream completed_file(completed_log_filename.c_str());
   
   // string line;
   // int line_number = 1;
   // while(getline(completed_file, line)) {
      
   //    if (line_number % 100000 == 0)
   //       cout << "On line number " << line_number << endl;

   //    istringstream iss(line);
   //    int event_number, run_number;
   //    iss >> run_number >> event_number;
      
   //    completed_events.insert(make_pair(to_string(run_number) + "_" + to_string(event_number), 1));

   //    line_number++;
   // }

   cout << "Finished loading \"completed\" events to memory." << endl;

   // ofstream completed_events_file_output;
   // completed_events_file_output.open(completed_log_filename.c_str(), ios::out | ios::app );


   vector<string> all_filenames;
   get_all_files_to_process(all_filenames, input_path);



   int file_counter = 0;
   // Loop through all those files and count events. 
   for (int i = 0; i < all_filenames.size(); i++) {

		ifstream file_to_process(all_filenames[i]);

		file_counter++;

		if ((file_counter % 100) == 0)
			cout << "Processing file number " << file_counter << " / " << all_filenames.size() << endl;

		MOD::Event event_being_read;
		int event_serial_number = 1;

		while (event_being_read.read_event(file_to_process)) {

			if ((event_serial_number % 10000) == 0)
			   cout << "Reading event number " << event_serial_number << endl;

			// count_events(event_being_read, trigger_numbers, trigger_prescales, npv, filtered_npv, jet_quality, passes_softdrop, assigned_trigger_fired, ak5_match, eta_24_cut);
			// cout << all_filenames[i] << endl;

			// Look up if we've already written this file.
			// auto search = completed_events.find( to_string(event_being_read.run_number()) + "_" + to_string(event_being_read.event_number()) );

			// if (search == completed_events.end()) {
				// Key not found. So proceed.

				// Find the name of the file.
				string file_for_current_event = find_correct_file_for_event(event_being_read.run_number(), event_being_read.event_number(), registry_info) + ".mod";

				// Open the file in append mode. 
				ofstream output_file;
   				output_file.open(output_path + "/" + file_for_current_event, ios::out | ios::app );

   				// Write the event to the file.
   				output_file << event_being_read;

				// Add it to the map completed_events and to the file.
	      		// completed_events.insert(make_pair(to_string(event_being_read.run_number()) + "_" + to_string(event_being_read.event_number()), 1));

				// completed_events_file_output << event_being_read.run_number() << "\t" << event_being_read.event_number() << endl; 
			}

			event_being_read = MOD::Event();
			event_serial_number++;
		// }
      
   }

  

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing everything in " << elapsed_seconds << " seconds!" << endl;

   return 0;
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





string find_correct_file_for_event(int run_number, int event_number, boost::unordered_map<string, string> & registry_info) {

   auto search = registry_info.find(to_string(run_number) + "_" + to_string(event_number));
   
   if(search != registry_info.end()) {
      // std::cout << "Found " << search->first << " " << search->second << '\n';
      return search->second;
   }
   else {
      //cout << "Oops didn't find it." << endl;
   }

   return "";
}

void load_registry_information(string registry_filename, boost::unordered_map<string, string> & registry_info) {

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