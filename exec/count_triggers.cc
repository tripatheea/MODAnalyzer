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
#include "fastjet/contrib/SoftDrop.hh"

#include "../interface/Event.h"

#include <boost/unordered_map.hpp>
#include <boost/filesystem.hpp>






// using namespace boost::filesystem;

using namespace std;
using namespace fastjet;
using namespace contrib;

using namespace boost::filesystem;


void count_triggers(MOD::Event & event_being_read, boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<string, int> & fired_trigger_numbers, boost::unordered_map<string, float> & fired_trigger_prescales);
void get_all_files_to_process(std::vector<string> & all_files, boost::filesystem::path input_path);

void write_stats(boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<string, int> & fired_trigger_numbers, boost::unordered_map<string, float> & fired_trigger_prescales);

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
   
   cout << endl << endl << "Starting count_triggers with the following given arguments: " << endl << endl;
   cout << "Input Path   : " << argv[1] << endl;
   cout << "Use LaTeX    : " << use_latex << endl;
   cout << "No. of Files : " << number_of_files_to_process << endl;
   cout << "No. of Events: " << number_of_events_to_process << endl;
   
   ofstream output_file("stats.dat", ios::out);

   // Recursively collect all filenames to process.

   vector<string> all_filenames;
   get_all_files_to_process(all_filenames, input_path);

   // Setup data structures to hold all counts. 
   boost::unordered_map<string, float> trigger_prescales;
   boost::unordered_map<string, int> trigger_numbers;
   boost::unordered_map<string, int> fired_trigger_numbers;
   boost::unordered_map<string, float> fired_trigger_prescales;
   
   int total_files = all_filenames.size();
   int total_validated_events = 0;

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

            if ((event_serial_number % 10000) == 0)
               cout << "Reading event number " << event_serial_number << endl;
            
            count_triggers(event_being_read, trigger_numbers, trigger_prescales, fired_trigger_numbers, fired_trigger_prescales);

            event_being_read = MOD::Event();
            event_serial_number++;
            total_validated_events++;
         }   
      }
      else {
         break;
      }

      write_stats(trigger_numbers, trigger_prescales, fired_trigger_numbers, fired_trigger_prescales);
      
   }

   
   cout << endl << "Everything done. Printing summary." << endl << endl;

   cout << endl << "==================================================================" << endl << endl;




   if (use_latex) {


      // cout << endl << endl << "================ LaTeX =============" << endl << endl;

   
      // cout << "\\hline" << endl;
      // cout << "\\hline" << endl;
      // cout << "& Events & Fraction\\\\" << endl;
      // cout << "\\hline" << endl;
      // cout << "Jet Primary Dataset & 20,022,826 & 1.000 \\\\" << endl;
      // cout << "Validated Run & " << total_validated_events << " & " << fixed << setprecision(3) << (total_validated_events/float(grand_total)) << " \\\\" << endl;
      // cout << "Loose Jet Quality (\\Tab{tab:jet_quality}) & " << jet_quality["LOOSE"] << " & " << fixed << setprecision(3) << (jet_quality["LOOSE"]/float(grand_total)) << " \\\\" << endl;
      // cout << "Assigned Trigger Fired (\\Tab{tab:trigger_table}) & " << assigned_trigger_fired << " & " << fixed << setprecision(3) << (assigned_trigger_fired/float(grand_total)) << " \\\\" << endl;
      // cout << "\\hline" << endl;
      // cout << "AK5 Match & " << ak5_match << " & " << fixed << setprecision(3) << (ak5_match/float(grand_total)) << "\\\\" << endl;
      // cout << "$\\left| \\eta \\right| < 2.4$ & " << eta_24_cut << " & " << fixed << setprecision(3) << (eta_24_cut/float(grand_total)) << " \\\\" << endl;
      // cout << "\\hline" << endl;
      // cout << "Passes Soft Drop ($z_g > z_{\\rm cut}$) & " << passes_softdrop[0.1] << " & " << fixed << setprecision(3) << (passes_softdrop[0.1]/float(grand_total)) << " \\\\" << endl;
      // cout << "\\hline" << endl;
      // cout << "\\hline" << endl;

      
      write_stats(trigger_numbers, trigger_prescales, fired_trigger_numbers, fired_trigger_prescales);
         

      
   }

   cout << endl << endl;
   

}


void write_stats(boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<string, int> & fired_trigger_numbers, boost::unordered_map<string, float> & fired_trigger_prescales) {

   ofstream output_file("./triggers_stats.dat", ios::out);

   for (auto kv : trigger_numbers) {
      output_file << kv.first << " \t  " << kv.second <<  " \t " << trigger_prescales[kv.first] << endl;
   }
   
   // // cout << endl << endl << "================ LaTeX =============" << endl << endl;

   // ofstream output_file("./stats.dat", ios::out);

   // // ostringstream LaTeX_stream;

   output_file << endl << endl << endl;

   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;
   output_file << "&Trigger & Trig Present & $\\langle$Prescale$\\rangle$ & Trig Fired? & $\\langle$Prescale$\\rangle$\\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "Single-jet & \\texttt{HLT\\_Jet15U} & " << trigger_numbers["HLT_Jet15U"] + trigger_numbers["HLT_Jet15U_v3"] << " & " << (trigger_prescales["HLT_Jet15U"] * trigger_numbers["HLT_Jet15U"] +  trigger_prescales["HLT_Jet15U_v3"] * trigger_numbers["HLT_Jet15U_v3"]) / (trigger_numbers["HLT_Jet15U"] + trigger_numbers["HLT_Jet15U_v3"]) << " & " << fired_trigger_numbers["HLT_Jet15U"] + fired_trigger_numbers["HLT_Jet15U_v3"] << " & " << (fired_trigger_prescales["HLT_Jet15U"] * fired_trigger_numbers["HLT_Jet15U"] +  fired_trigger_prescales["HLT_Jet15U_v3"] * fired_trigger_numbers["HLT_Jet15U_v3"]) / (fired_trigger_numbers["HLT_Jet15U"] + fired_trigger_numbers["HLT_Jet15U_v3"]) << " \\\\" << endl;
   output_file << "&\\texttt{* HLT\\_Jet15U\\_HNF} & " << trigger_numbers["HLT_Jet15U_HcalNoiseFiltered"] + trigger_numbers["HLT_Jet15U_HcalNoiseFiltered_v3"] << " & " << (trigger_prescales["HLT_Jet15U_HcalNoiseFiltered"] * trigger_numbers["HLT_Jet15U_HcalNoiseFiltered"] +  trigger_prescales["HLT_Jet15U_HcalNoiseFiltered_v3"] * trigger_numbers["HLT_Jet15U_HcalNoiseFiltered_v3"]) / (trigger_numbers["HLT_Jet15U_HcalNoiseFiltered"] + trigger_numbers["HLT_Jet15U_HcalNoiseFiltered_v3"]) << " & " << fired_trigger_numbers["HLT_Jet15U_HcalNoiseFiltered"] + fired_trigger_numbers["HLT_Jet15U_HcalNoiseFiltered_v3"] << " & " << (fired_trigger_prescales["HLT_Jet15U_HcalNoiseFiltered"] * fired_trigger_numbers["HLT_Jet15U_HcalNoiseFiltered"] +  fired_trigger_prescales["HLT_Jet15U_HcalNoiseFiltered_v3"] * fired_trigger_numbers["HLT_Jet15U_HcalNoiseFiltered_v3"]) / (fired_trigger_numbers["HLT_Jet15U_HcalNoiseFiltered"] + fired_trigger_numbers["HLT_Jet15U_HcalNoiseFiltered_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet30U} & " << trigger_numbers["HLT_Jet30U"] + trigger_numbers["HLT_Jet30U_v3"] << " & " << (trigger_prescales["HLT_Jet30U"] * trigger_numbers["HLT_Jet30U"] +  trigger_prescales["HLT_Jet30U_v3"] * trigger_numbers["HLT_Jet30U_v3"]) / (trigger_numbers["HLT_Jet30U"] + trigger_numbers["HLT_Jet30U_v3"]) << " & " << fired_trigger_numbers["HLT_Jet30U"] + fired_trigger_numbers["HLT_Jet30U_v3"] << " & " << (fired_trigger_prescales["HLT_Jet30U"] * fired_trigger_numbers["HLT_Jet30U"] +  fired_trigger_prescales["HLT_Jet30U_v3"] * fired_trigger_numbers["HLT_Jet30U_v3"]) / (fired_trigger_numbers["HLT_Jet30U"] + fired_trigger_numbers["HLT_Jet30U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet50U} & " << trigger_numbers["HLT_Jet50U"] + trigger_numbers["HLT_Jet50U_v3"] << " & " << (trigger_prescales["HLT_Jet50U"] * trigger_numbers["HLT_Jet50U"] +  trigger_prescales["HLT_Jet50U_v3"] * trigger_numbers["HLT_Jet50U_v3"]) / (trigger_numbers["HLT_Jet50U"] + trigger_numbers["HLT_Jet50U_v3"]) << " &  " << fired_trigger_numbers["HLT_Jet50U"] + fired_trigger_numbers["HLT_Jet50U_v3"] << " & " << (fired_trigger_prescales["HLT_Jet50U"] * fired_trigger_numbers["HLT_Jet50U"] +  fired_trigger_prescales["HLT_Jet50U_v3"] * fired_trigger_numbers["HLT_Jet50U_v3"]) / (fired_trigger_numbers["HLT_Jet50U"] + fired_trigger_numbers["HLT_Jet50U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet70U} & " << trigger_numbers["HLT_Jet70U"] + trigger_numbers["HLT_Jet70U_v2"] + trigger_numbers["HLT_Jet70U_v3"] << " & " << (trigger_prescales["HLT_Jet70U"] * trigger_numbers["HLT_Jet70U"] + trigger_prescales["HLT_Jet70U_v2"] * trigger_numbers["HLT_Jet70U_v2"] +  trigger_prescales["HLT_Jet70U_v3"] * trigger_numbers["HLT_Jet70U_v3"]) / (trigger_numbers["HLT_Jet70U"] + trigger_numbers["HLT_Jet70U_v2"] + trigger_numbers["HLT_Jet70U_v3"]) << " & " << fired_trigger_numbers["HLT_Jet70U"] + fired_trigger_numbers["HLT_Jet70U_v2"] + fired_trigger_numbers["HLT_Jet70U_v3"] << " & " << (fired_trigger_prescales["HLT_Jet70U"] * fired_trigger_numbers["HLT_Jet70U"] + fired_trigger_prescales["HLT_Jet70U_v2"] * fired_trigger_numbers["HLT_Jet70U_v2"] +  fired_trigger_prescales["HLT_Jet70U_v3"] * fired_trigger_numbers["HLT_Jet70U_v3"]) / (fired_trigger_numbers["HLT_Jet70U"] + fired_trigger_numbers["HLT_Jet70U_v2"] + fired_trigger_numbers["HLT_Jet70U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet100U} & " << trigger_numbers["HLT_Jet100U"] + trigger_numbers["HLT_Jet100U_v2"] + trigger_numbers["HLT_Jet100U_v3"] << " & " << (trigger_prescales["HLT_Jet100U"] * trigger_numbers["HLT_Jet100U"] + trigger_prescales["HLT_Jet100U_v2"] * trigger_numbers["HLT_Jet100U_v2"] +  trigger_prescales["HLT_Jet100U_v3"] * trigger_numbers["HLT_Jet100U_v3"]) / (trigger_numbers["HLT_Jet100U"] + trigger_numbers["HLT_Jet100U_v2"] + trigger_numbers["HLT_Jet100U_v3"]) << "  & " << fired_trigger_numbers["HLT_Jet100U"] + fired_trigger_numbers["HLT_Jet100U_v2"] + fired_trigger_numbers["HLT_Jet100U_v3"] << " & " << (fired_trigger_prescales["HLT_Jet100U"] * fired_trigger_numbers["HLT_Jet100U"] + fired_trigger_prescales["HLT_Jet100U_v2"] * fired_trigger_numbers["HLT_Jet100U_v2"] +  fired_trigger_prescales["HLT_Jet100U_v3"] * fired_trigger_numbers["HLT_Jet100U_v3"]) / (fired_trigger_numbers["HLT_Jet100U"] + fired_trigger_numbers["HLT_Jet100U_v2"] + fired_trigger_numbers["HLT_Jet100U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet140U} & " << trigger_numbers["HLT_Jet140U_v1"] + trigger_numbers["HLT_Jet140U_v3"] << " & " << (trigger_prescales["HLT_Jet140U_v1"] * trigger_numbers["HLT_Jet140U_v1"] +  trigger_prescales["HLT_Jet140U_v3"] * trigger_numbers["HLT_Jet140U_v3"]) / (trigger_numbers["HLT_Jet140U_v1"] + trigger_numbers["HLT_Jet140U_v3"]) << "  & " << fired_trigger_numbers["HLT_Jet140U_v1"] + fired_trigger_numbers["HLT_Jet140U_v3"] << " & " << (fired_trigger_prescales["HLT_Jet140U_v1"] * fired_trigger_numbers["HLT_Jet140U_v1"] +  fired_trigger_prescales["HLT_Jet140U_v3"] * fired_trigger_numbers["HLT_Jet140U_v3"]) / (fired_trigger_numbers["HLT_Jet140U_v1"] + fired_trigger_numbers["HLT_Jet140U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet180U} & " << trigger_numbers["HLT_Jet180U"] << " & " << trigger_prescales["HLT_Jet180U"] << "  & " << fired_trigger_numbers["HLT_Jet180U"] << " & " << fired_trigger_prescales["HLT_Jet180U"] << "  \\\\" << endl;   
   
   output_file << "\\hline" << endl;

   output_file << "Di-jet & \\texttt{HLT\\_DiJetAve15U} & " << trigger_numbers["HLT_DiJetAve15U"] + trigger_numbers["HLT_DiJetAve15U_v3"] << " & " << (trigger_prescales["HLT_DiJetAve15U"] * trigger_numbers["HLT_DiJetAve15U"] +  trigger_prescales["HLT_DiJetAve15U_v3"] * trigger_numbers["HLT_DiJetAve15U_v3"]) / (trigger_numbers["HLT_DiJetAve15U"] + trigger_numbers["HLT_DiJetAve15U_v3"]) << " \\\\" << endl;
   output_file << "&\\texttt{* HLT\\_DiJetAve30U} & " << trigger_numbers["HLT_DiJetAve30U"] + trigger_numbers["HLT_DiJetAve30U_v3"] << " & " << (trigger_prescales["HLT_DiJetAve30U"] * trigger_numbers["HLT_DiJetAve30U"] +  trigger_prescales["HLT_DiJetAve30U_v3"] * trigger_numbers["HLT_DiJetAve30U_v3"]) / (trigger_numbers["HLT_DiJetAve30U"] + trigger_numbers["HLT_DiJetAve30U_v3"]) << "  & " << fired_trigger_numbers["HLT_DiJetAve30U"] + fired_trigger_numbers["HLT_DiJetAve30U_v3"] << " & " << (fired_trigger_prescales["HLT_DiJetAve30U"] * fired_trigger_numbers["HLT_DiJetAve30U"] +  fired_trigger_prescales["HLT_DiJetAve30U_v3"] * fired_trigger_numbers["HLT_DiJetAve30U_v3"]) / (fired_trigger_numbers["HLT_DiJetAve30U"] + fired_trigger_numbers["HLT_DiJetAve30U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_DiJetAve50U} & " << trigger_numbers["HLT_DiJetAve50U"] + trigger_numbers["HLT_DiJetAve50U_v3"] << " & " << (trigger_prescales["HLT_DiJetAve50U"] * trigger_numbers["HLT_DiJetAve50U"] +  trigger_prescales["HLT_DiJetAve50U_v3"] * trigger_numbers["HLT_DiJetAve50U_v3"]) / (trigger_numbers["HLT_DiJetAve50U"] + trigger_numbers["HLT_DiJetAve50U_v3"]) << " & " << fired_trigger_numbers["HLT_DiJetAve50U"] + fired_trigger_numbers["HLT_DiJetAve50U_v3"] << " & " << (fired_trigger_prescales["HLT_DiJetAve50U"] * fired_trigger_numbers["HLT_DiJetAve50U"] +  fired_trigger_prescales["HLT_DiJetAve50U_v3"] * fired_trigger_numbers["HLT_DiJetAve50U_v3"]) / (fired_trigger_numbers["HLT_DiJetAve50U"] + fired_trigger_numbers["HLT_DiJetAve50U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_DiJetAve70U} & " << trigger_numbers["HLT_DiJetAve70U"] + trigger_numbers["HLT_DiJetAve70U_v2"] + trigger_numbers["HLT_DiJetAve70U_v3"] << " & " << (trigger_prescales["HLT_DiJetAve70U"] * trigger_numbers["HLT_DiJetAve70U"] + trigger_prescales["HLT_DiJetAve70U_v2"] * trigger_numbers["HLT_DiJetAve70U_v2"] +  trigger_prescales["HLT_DiJetAve70U_v3"] * trigger_numbers["HLT_DiJetAve70U_v3"]) / (trigger_numbers["HLT_DiJetAve70U"] + trigger_numbers["HLT_DiJetAve70U_v2"] + trigger_numbers["HLT_DiJetAve70U_v3"]) << "  & " << fired_trigger_numbers["HLT_DiJetAve70U"] + fired_trigger_numbers["HLT_DiJetAve70U_v2"] + fired_trigger_numbers["HLT_DiJetAve70U_v3"] << " & " << (fired_trigger_prescales["HLT_DiJetAve70U"] * fired_trigger_numbers["HLT_DiJetAve70U"] + fired_trigger_prescales["HLT_DiJetAve70U_v2"] * fired_trigger_numbers["HLT_DiJetAve70U_v2"] +  fired_trigger_prescales["HLT_DiJetAve70U_v3"] * fired_trigger_numbers["HLT_DiJetAve70U_v3"]) / (fired_trigger_numbers["HLT_DiJetAve70U"] + fired_trigger_numbers["HLT_DiJetAve70U_v2"] + fired_trigger_numbers["HLT_DiJetAve70U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_DiJetAve100U} & " << trigger_numbers["HLT_DiJetAve100U_v1"] + trigger_numbers["HLT_DiJetAve100U_v3"] << " & " << (trigger_prescales["HLT_DiJetAve100U_v1"] * trigger_numbers["HLT_DiJetAve100U_v1"] + trigger_prescales["HLT_DiJetAve100U_v2"] * trigger_numbers["HLT_DiJetAve100U_v2"] +  trigger_prescales["HLT_DiJetAve100U_v3"] * trigger_numbers["HLT_DiJetAve100U_v3"]) / (trigger_numbers["HLT_DiJetAve100U"] + trigger_numbers["HLT_DiJetAve100U_v2"] + trigger_numbers["HLT_DiJetAve100U_v3"]) << "  & " << fired_trigger_numbers["HLT_DiJetAve100U_v1"] + fired_trigger_numbers["HLT_DiJetAve100U_v3"] << " & " << (fired_trigger_prescales["HLT_DiJetAve100U_v1"] * fired_trigger_numbers["HLT_DiJetAve100U_v1"] + fired_trigger_prescales["HLT_DiJetAve100U_v2"] * fired_trigger_numbers["HLT_DiJetAve100U_v2"] +  fired_trigger_prescales["HLT_DiJetAve100U_v3"] * fired_trigger_numbers["HLT_DiJetAve100U_v3"]) / (fired_trigger_numbers["HLT_DiJetAve100U"] + fired_trigger_numbers["HLT_DiJetAve100U_v2"] + fired_trigger_numbers["HLT_DiJetAve100U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_DiJetAve100U} & " << trigger_numbers["HLT_DiJetAve100U_v1"] + trigger_numbers["HLT_DiJetAve100U_v3"] << " & " << (trigger_prescales["HLT_DiJetAve100U_v1"] * trigger_numbers["HLT_DiJetAve100U_v1"] + trigger_prescales["HLT_DiJetAve100U_v2"] * trigger_numbers["HLT_DiJetAve100U_v2"] +  trigger_prescales["HLT_DiJetAve100U_v3"] * trigger_numbers["HLT_DiJetAve100U_v3"]) / (trigger_numbers["HLT_DiJetAve100U"] + trigger_numbers["HLT_DiJetAve100U_v2"] + trigger_numbers["HLT_DiJetAve100U_v3"]) << "  & " << fired_trigger_numbers["HLT_DiJetAve100U_v1"] + fired_trigger_numbers["HLT_DiJetAve100U_v3"] << " & " << (fired_trigger_prescales["HLT_DiJetAve100U_v1"] * fired_trigger_numbers["HLT_DiJetAve100U_v1"] + fired_trigger_prescales["HLT_DiJetAve100U_v2"] * fired_trigger_numbers["HLT_DiJetAve100U_v2"] +  fired_trigger_prescales["HLT_DiJetAve100U_v3"] * fired_trigger_numbers["HLT_DiJetAve100U_v3"]) / (fired_trigger_numbers["HLT_DiJetAve100U"] + fired_trigger_numbers["HLT_DiJetAve100U_v2"] + fired_trigger_numbers["HLT_DiJetAve100U_v3"]) << "  \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_DiJetAve140U} & " << trigger_numbers["HLT_DiJetAve140U"] << " & " << trigger_prescales["HLT_DiJetAve140U"] <<  " & " << fired_trigger_numbers["HLT_DiJetAve140U"] << " & " << fired_trigger_prescales["HLT_DiJetAve140U"] << "  \\\\" << endl;   
   
   output_file << "\\hline" << endl;

   output_file << "Quad-jet & \\texttt{HLT\\_QuadJet20U} & " << trigger_numbers["HLT_QuadJet20U"] << " & " << trigger_prescales["HLT_QuadJet20U"] << " & " << fired_trigger_numbers["HLT_QuadJet20U"] << " & " << fired_trigger_prescales["HLT_QuadJet20U"] << " \\\\" << endl;
   output_file << "& \\texttt{HLT\\_QuadJet25U} & " << trigger_numbers["HLT_QuadJet25U"] << " & " << trigger_prescales["HLT_QuadJet25U"] << " & " << fired_trigger_numbers["HLT_QuadJet25U"] << " & " << fired_trigger_prescales["HLT_QuadJet25U"] << " \\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "$H_T$ & \\texttt{HLT\\_HT100U} & " << trigger_numbers["HLT_HT100U"] << " & " << trigger_prescales["HLT_HT100U"] << " & " << fired_trigger_numbers["HLT_HT100U"] << " & " << fired_trigger_prescales["HLT_HT100U"] << " \\\\" << endl;
   output_file << "&\\texttt{HLT\\_HT120U} & " << trigger_numbers["HLT_HT120U"] <<  " & " << trigger_prescales["HLT_HT120U"] <<  " & " << fired_trigger_numbers["HLT_HT120U"] <<  " & " << fired_trigger_prescales["HLT_HT120U"] << " \\\\" << endl;
   output_file << "&\\texttt{HLT\\_HT140U} & " << trigger_numbers["HLT_HT140U"] <<  " & " << trigger_prescales["HLT_HT140U"] << " & " << fired_trigger_numbers["HLT_HT140U"] <<  " & " << fired_trigger_prescales["HLT_HT140U"] << " \\\\" << endl;
   output_file << "&\\texttt{HLT\\_EcalOnly\\_SumEt160} & " << trigger_numbers["HLT_EcalOnly_SumEt160"] << " & " << trigger_prescales["HLT_EcalOnly_SumEt160"] << " & " << fired_trigger_numbers["HLT_EcalOnly_SumEt160"] << " & " << fired_trigger_prescales["HLT_EcalOnly_SumEt160"] << " \\\\" << endl;

   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;


   output_file << endl << endl << endl << endl << endl;

   // for (auto kv : npv) {
   //    output_file << kv.first << " => " << kv.second << endl;
   // }
   // output_file << endl << endl;

   // // cout << LaTeX_stream.str() << endl;   

   // // output_file << LaTeX_stream;   
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

void count_triggers(MOD::Event & event_being_read, boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<string, int> & fired_trigger_numbers, boost::unordered_map<string, float> & fired_trigger_prescales) {
   
   vector<MOD::Trigger> all_triggers = event_being_read.triggers();

   for (unsigned i=0; i < all_triggers.size(); i++) {
      MOD::Trigger trigger = all_triggers[i];


      // First present triggers.
      auto numbers_search = trigger_numbers.find(trigger.name());

      if (numbers_search != trigger_numbers.end()) {
         numbers_search->second++;
      }
      else {
         trigger_numbers.insert(make_pair(trigger.name(), 1));
      }

      // Average prescale.
      auto prescale_search = trigger_prescales.find(trigger.name());

      if (prescale_search != trigger_prescales.end()) {
         // Get the total number of prescales already in the hashmap.
         int n = trigger_numbers[trigger.name()];
         float summation_x = n * trigger_prescales[trigger.name()];
         float new_mean = (summation_x + trigger.prescale()) / (n + 1);

         trigger_prescales[trigger.name()] = new_mean;   
      }
      else {
         trigger_prescales.insert(make_pair(trigger.name(), trigger.prescale()));
      }


      // Next, fired triggers.
      if (trigger.fired()) {
         auto numbers_search = fired_trigger_numbers.find(trigger.name());

         if (numbers_search != fired_trigger_numbers.end()) {
            numbers_search->second++;
         }
         else {
            fired_trigger_numbers.insert(make_pair(trigger.name(), 1));
         }

         // Average prescale.
         auto prescale_search = fired_trigger_prescales.find(trigger.name());

         if (prescale_search != fired_trigger_prescales.end()) {
            // Get the total number of prescales already in the hashmap.
            int n = fired_trigger_numbers[trigger.name()];
            float summation_x = n * fired_trigger_prescales[trigger.name()];
            float new_mean = (summation_x + trigger.prescale()) / (n + 1);

            fired_trigger_prescales[trigger.name()] = new_mean;   
         }
         else {
            fired_trigger_prescales.insert(make_pair(trigger.name(), trigger.prescale()));
         }
      }
      
   }


   // PseudoJet trigger_jet = event_being_read.trigger_jet();
   // event_being_read.convert_to_one_jet();


   // if (event_being_read.trigger_jet_is_matched() && (trigger_jet.user_info<MOD::InfoCalibratedJet>().jet_quality() >= 1)) {   // Jet quality level: FAILED = 0, LOOSE = 1, MEDIUM = 2, TIGHT = 3      
   //    output_file << event_being_read;
   // }
   
}
