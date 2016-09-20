#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>


#include <boost/unordered_map.hpp>
#include <chrono>


using namespace std;


int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();



   string input_filename(argv[1]);


   ifstream input_file(input_filename.c_str());
   
   cout << endl << endl << "Starting process with the following given arguments: " << endl;
   cout << "Input File: " << argv[1] << endl << endl;

   int line_number = 0;
   int duplicates = 0;
   int unique = 0;

   boost::unordered_map<pair<int, int>, int> all_events;

   string line;
   // while ((getline(input_file, line)) && (line_number < 10000)) {
   while (getline(input_file, line)) {

      if (line_number % 100000 == 0)
         cout << "On line number " << line_number << endl;

      istringstream iss(line);

      int event_number, run_number;
      
      iss >> run_number >> event_number;

      if (all_events.find(make_pair(run_number, event_number)) == all_events.end()) {
         all_events.insert(make_pair(make_pair(run_number, event_number), 1));
         unique++;
      }
      else {
         // Found a duplicate.
         duplicates++;
      }
      

      line_number++;
   }

   cout << endl;
   cout << "Parsed " << line_number << " events." << endl;
   cout << "Found " << duplicates << " duplicate events." << endl;
   cout << "Found " << unique << " unique events." << endl << endl;

   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << line_number << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}

