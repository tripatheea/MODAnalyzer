#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <iterator>
#include <iomanip>
#include <limits>
#include <chrono>
#include <boost/unordered_map.hpp>

using namespace std;


int main(int argc, char * argv[]) {

   auto start = std::chrono::steady_clock::now();

   if (argc <= 2) {
        std::cerr << "ERROR: You need to supply five arguments- first, path to the input data; second, path to the output file; third, number of events to process. The path has to be either absolute or relative to the bin directory:" << std::endl << std::endl << "./analysis (input_file.dat) (output_file.dat) [optional Nev]" << std::endl;
        return 1;
   }

   ifstream registry_file(argv[1]);
   ofstream output_file(argv[2], ios::out);


   cout << endl << endl << "Starting analysis with the following given arguments: " << endl;
   cout << "Registry file: " << argv[1] << endl;
   cout << "Output file: " << argv[2] << endl;

   

   boost::unordered_map<string, int> file_numbers;

   int event_serial_number = 1;
   string line;
   while(getline(registry_file, line)) {
      string path, filename;
      int event_number, run_number;

      istringstream iss(line);

      iss >> event_number >> run_number >> path >> filename;

      filename = filename.substr(0, filename.length() - 5);

      auto search = file_numbers.find(filename);
      if (search != file_numbers.end()) {
         search->second++;
      }
      else {
         file_numbers.insert(make_pair(filename, 1));
      }

      event_serial_number++;
   }

   vector<string> filenames;
   for(auto kv : file_numbers) {
      filenames.push_back(kv.first);
   }

   sort( filenames.begin(), filenames.end() );


   for (string file : filenames) {
      output_file << file << "\t" << file_numbers[file] << endl;
   }


   auto finish = std::chrono::steady_clock::now();
   double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start).count();
   cout << "Finished processing " << (event_serial_number - 1) << " events in " << elapsed_seconds << " seconds!" << endl;

   return 0;
}




