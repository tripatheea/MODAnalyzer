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
#include <boost/regex.hpp>


// test



// using namespace boost::filesystem;

using namespace std;
using namespace fastjet;
using namespace contrib;

using namespace boost::filesystem;


void count_events(MOD::Event & event_being_read, boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<int, int> & npv, boost::unordered_map<int, int> & filtered_npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut);
void get_all_files_to_process(std::vector<string> & all_files, vector<boost::filesystem::path> input_paths);
void count_triggers(MOD::Event & event_being_read, boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<string, int> & fired_trigger_numbers, boost::unordered_map<string, float> & fired_trigger_prescales);



void write_stats(string filename, boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<int, int> & npv, boost::unordered_map<int, int> & filtered_npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut, int & total_validated_events, int & grand_total, boost::unordered_map<string, int> & all_trigger_numbers, boost::unordered_map<string, float> & all_trigger_prescales, boost::unordered_map<string, int> & fired_all_trigger_numbers, boost::unordered_map<string, float> & fired_all_trigger_prescales);

int main(int argc, char * argv[]) {
   
   
   string filename(argv[1]);
   
   
   if (argc < 3) {
      std::cerr << "You need to give at least two arguments: output path and one or more input directories." << endl;
      return 0;
   }

   vector<path> input_paths;
   for (unsigned i=2; i < argc; i++) {
     input_paths.push_back(argv[i]);
   }


   
   cout << endl << endl << "Starting count_events with the following given arguments: " << endl << endl;
   cout << "Output File     : " << argv[1] << endl;
   cout << "Input Paths (" << input_paths.size() << ") : ";
   for (unsigned i=0; i < input_paths.size(); i++) {
    if (i > 0)
        cout << "                  " << input_paths[i] << endl;
    else 
        cout << input_paths[i] << endl;
   }   

   // Recursively collect all filenames to process.

   vector<string> all_filenames;
   get_all_files_to_process(all_filenames, input_paths);

   // Sort the list.
   std::sort(all_filenames.begin(), all_filenames.end());
   
   cout << endl << endl << "Starting counting event on " << all_filenames.size() << " files." << endl << endl;

   // Setup data structures to hold all counts. 
   boost::unordered_map<int, int> npv;
   boost::unordered_map<int, int> filtered_npv;
   boost::unordered_map<string, int> jet_quality;
   boost::unordered_map<float, int> passes_softdrop;

   boost::unordered_map<string, int> trigger_numbers;
   boost::unordered_map<string, float> trigger_prescales;


   // Setup data structures to hold all trigger counts. 
   boost::unordered_map<string, float> all_trigger_prescales;
   boost::unordered_map<string, int> all_trigger_numbers;
   boost::unordered_map<string, int> fired_all_trigger_numbers;
   boost::unordered_map<string, float> fired_all_trigger_prescales;



   // int total_files = all_filenames.size();
   int total_files = all_filenames.size();
   int total_validated_events = 0;
   int assigned_trigger_fired = 0;
   int ak5_match = 0;
   int eta_24_cut = 0;


   // Remove the following line.
   ofstream output_file(filename, ios::out);
   

   int grand_total = 20022826;
   // int grand_total = 20022;

   int file_counter = 0;
   // Loop through all those files and count events. 


   for (int i = 0; i < total_files; i++) {

      
         ifstream file_to_process(all_filenames[i]);

         // cout << endl << (file_counter + 1) << "/" << total_files << "\t" << all_filenames[i] << endl;
         
         file_counter++;

         if ((file_counter % 100) == 0)
            cout << "Processing file number " << file_counter << " / " << total_files << endl;

         MOD::Event event_being_read;
         int event_serial_number = 1;

         // cout << all_filenames[i] << endl;

         while (event_being_read.read_event(file_to_process)) {

            // if ((event_serial_number % 10000) == 0)
               // cout << "Reading event number " << event_serial_number << endl;
            
            count_events(event_being_read, trigger_numbers, trigger_prescales, npv, filtered_npv, jet_quality, passes_softdrop, assigned_trigger_fired, ak5_match, eta_24_cut);
            count_triggers(event_being_read, all_trigger_numbers, all_trigger_prescales, fired_all_trigger_numbers, fired_all_trigger_prescales);

            event_being_read = MOD::Event();
            event_serial_number++;
            total_validated_events++;
         }

         
        

         cout << "Total Validated Events: " << total_validated_events << endl;   
        

      

      write_stats(filename, trigger_numbers, trigger_prescales, npv, filtered_npv, jet_quality, passes_softdrop, assigned_trigger_fired, ak5_match, eta_24_cut, total_validated_events, grand_total, all_trigger_numbers, all_trigger_prescales, fired_all_trigger_numbers, fired_all_trigger_prescales);
      
   }

   
   cout << endl << "Everything done. Printing summary." << endl << endl;

   cout << endl << "==================================================================" << endl << endl;


    write_stats(filename, trigger_numbers, trigger_prescales, npv, filtered_npv, jet_quality, passes_softdrop, assigned_trigger_fired, ak5_match, eta_24_cut, total_validated_events, grand_total, all_trigger_numbers, all_trigger_prescales, fired_all_trigger_numbers, fired_all_trigger_prescales);
         

      

   

}













void write_stats(string filename, boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<int, int> & npv, boost::unordered_map<int, int> & filtered_npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut, int & total_validated_events, int & grand_total, boost::unordered_map<string, int> & all_trigger_numbers, boost::unordered_map<string, float> & all_trigger_prescales, boost::unordered_map<string, int> & fired_all_trigger_numbers, boost::unordered_map<string, float> & fired_all_trigger_prescales) {
   // cout << endl << endl << "================ LaTeX =============" << endl << endl;

   ofstream output_file(filename, ios::out);

   

   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;
   output_file << "& Events & Fraction\\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "Jet Primary Dataset & 20,022,826 & 1.000 \\\\" << endl;
   output_file << "Validated Run & " << total_validated_events << " & " << fixed << setprecision(3) << (total_validated_events/float(grand_total)) << " \\\\" << endl;
   output_file << "Assigned Trigger Fired (\\Tab{tab:trigger_table}) & " << assigned_trigger_fired << " & " << fixed << setprecision(3) << (assigned_trigger_fired/float(grand_total)) << " \\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "Loose Jet Quality (\\Tab{tab:jet_quality}) & " << jet_quality["LOOSE"] << " & " << fixed << setprecision(3) << (jet_quality["LOOSE"]/float(grand_total)) << " \\\\" << endl;
   output_file << "AK5 Match & " << ak5_match << " & " << fixed << setprecision(3) << (ak5_match/float(grand_total)) << "\\\\" << endl;
   output_file << "$\\left| \\eta \\right| < 2.4$ & " << eta_24_cut << " & " << fixed << setprecision(3) << (eta_24_cut/float(grand_total)) << " \\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "Passes Soft Drop ($z_g > z_{\\rm cut}$) & " << passes_softdrop[0.1] << " & " << fixed << setprecision(3) << (passes_softdrop[0.1]/float(grand_total)) << " \\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;


   output_file << endl << endl << endl << endl << endl;

   float npv_total = 0;
   float npv_larger_than_5_total = 0;
   for (auto kv : npv) {
      // output_file << kv.first << " => " << kv.second << endl;
      npv_total += kv.second;
      if (kv.first >= 15)
         npv_larger_than_5_total += kv.second;
   }

   float filtered_npv_total = 0;
   float filtered_npv_larger_than_5_total = 0;
   for (auto kv : filtered_npv) {
      // output_file << kv.first << " => " << kv.second << endl;
      filtered_npv_total += kv.second;
      if (kv.first >= 15)
         filtered_npv_larger_than_5_total += kv.second;
   }


   output_file << endl << endl;

   

   // NPV.
   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;
   output_file << "& All & & Filtered \\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "$N_{\\rm PV}$ & Events & Fraction & Events & Fraction\\\\" << endl;
   output_file << "\\hline" << endl;\

   output_file << "1 & " << npv[1] << " & " << (npv[1] / npv_total) << " & " << filtered_npv[1] << " & " << (filtered_npv[1] / filtered_npv_total) << "\\\\" << endl;
   output_file << "2 & " << npv[2] << " & " << (npv[2] / npv_total) << " & " << filtered_npv[2] << " & " << (filtered_npv[2] / filtered_npv_total) << "\\\\" << endl;
   output_file << "3 & " << npv[3] << " & " << (npv[3] / npv_total) << " & " << filtered_npv[3] << " & " << (filtered_npv[3] / filtered_npv_total) << "\\\\" << endl;
   output_file << "4 & " << npv[4] << " & " << (npv[4] / npv_total) << " & " << filtered_npv[4] << " & " << (filtered_npv[4] / filtered_npv_total) << "\\\\" << endl;
   output_file << "5 & " << npv[5] << " & " << (npv[5] / npv_total) << " & " << filtered_npv[5] << " & " << (filtered_npv[5] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "6 & " << npv[6] << " & " << (npv[6] / npv_total) << " & " << filtered_npv[6] << " & " << (filtered_npv[6] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "7 & " << npv[7] << " & " << (npv[7] / npv_total) << " & " << filtered_npv[7] << " & " << (filtered_npv[7] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "8 & " << npv[8] << " & " << (npv[8] / npv_total) << " & " << filtered_npv[8] << " & " << (filtered_npv[8] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "9 & " << npv[9] << " & " << (npv[9] / npv_total) << " & " << filtered_npv[9] << " & " << (filtered_npv[9] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "10 & " << npv[10] << " & " << (npv[10] / npv_total) << " & " << filtered_npv[10] << " & " << (filtered_npv[10] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "11 & " << npv[11] << " & " << (npv[11] / npv_total) << " & " << filtered_npv[11] << " & " << (filtered_npv[11] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "12 & " << npv[12] << " & " << (npv[12] / npv_total) << " & " << filtered_npv[12] << " & " << (filtered_npv[12] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "13 & " << npv[13] << " & " << (npv[13] / npv_total) << " & " << filtered_npv[13] << " & " << (filtered_npv[13] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "14 & " << npv[14] << " & " << (npv[14] / npv_total) << " & " << filtered_npv[14] << " & " << (filtered_npv[14] / filtered_npv_total) << "\\\\" << endl;   
   output_file << "$\\ge 15$ & " << npv_larger_than_5_total << " & " << (npv_larger_than_5_total / npv_total) << " & " << filtered_npv_larger_than_5_total << " & " << (filtered_npv_larger_than_5_total / filtered_npv_total) << "\\\\" << endl;   
   
   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;

   // cout << LaTeX_stream.str() << endl;   

   output_file << endl << endl << endl << endl << endl;

   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;
   output_file << "Hardest Jet $p_T$ & Trigger Name & Events  & $\\langle$Prescale$\\rangle$ \\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "$[85, 115]~\\GeV$ & \\texttt{HLT\\_Jet30U} & " << trigger_numbers["Jet30U"] <<  " & " << trigger_prescales["Jet30U"]  << " \\\\" << endl;
   output_file << "$[115, 150]~\\GeV$ & \\texttt{HLT\\_Jet50U} & " << trigger_numbers["Jet50U"] << " & " << trigger_prescales["Jet50U"] << " \\\\" << endl;
   output_file << "$[150, 200]~\\GeV$ & \\texttt{HLT\\_Jet70U}  & " << trigger_numbers["Jet70U"] << " & " << trigger_prescales["Jet70U"] << " \\\\" << endl;
   output_file << "$[200, 250]~\\GeV$ & \\texttt{HLT\\_Jet100U} & " << trigger_numbers["Jet100U"] << " & " << trigger_prescales["Jet100U"] <<  "\\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "\\multirow{2}{*}{$> 250~\\GeV$} & \\texttt{HLT\\_Jet100U}   & " << trigger_numbers["extra_100U"] << " & " << trigger_prescales["extra_100U"] << "\\\\" << endl;
   output_file << "& \\texttt{HLT\\_Jet140U} & " << trigger_numbers["Jet140U"] << " & " << trigger_prescales["Jet140U"] << "\\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;

   // output_file << LaTeX_stream;      


   output_file << endl << endl << endl << endl << endl;



   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;
   output_file << "&Trigger & Trig Present & $\\langle$Prescale$\\rangle$ & Trig Fired? & $\\langle$Prescale$\\rangle$\\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "Single-jet & \\texttt{HLT\\_Jet15U} & "   << all_trigger_numbers["Jet15U"]                   << " & " << all_trigger_prescales["Jet15U"]                   << " & " << fired_all_trigger_numbers["Jet15U"]                   << " & " << fired_all_trigger_prescales["Jet15U"]                   << " \\\\" << endl;
   output_file << "&\\texttt{* HLT\\_Jet15U\\_HNF} & "       << all_trigger_numbers["Jet15U_HcalNoiseFiltered"] << " & " << all_trigger_prescales["Jet15U_HcalNoiseFiltered"] << " & " << fired_all_trigger_numbers["Jet15U_HcalNoiseFiltered"] << " & " << fired_all_trigger_prescales["Jet15U_HcalNoiseFiltered"] << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet30U} & "             << all_trigger_numbers["Jet30U"]                   << " & " << all_trigger_prescales["Jet30U"]                   << " & " << fired_all_trigger_numbers["Jet30U"]                   << " & " << fired_all_trigger_prescales["Jet30U"]                   << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet50U} & "             << all_trigger_numbers["Jet50U"]                   << " & " << all_trigger_prescales["Jet50U"]                   << " & " << fired_all_trigger_numbers["Jet50U"]                   << " & " << fired_all_trigger_prescales["Jet50U"]                   << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet70U} & "             << all_trigger_numbers["Jet70U"]                   << " & " << all_trigger_prescales["Jet70U"]                   << " & " << fired_all_trigger_numbers["Jet70U"]                   << " & " << fired_all_trigger_prescales["Jet70U"]                   << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet100U} & "            << all_trigger_numbers["Jet100U"]                  << " & " << all_trigger_prescales["Jet100U"]                  << " & " << fired_all_trigger_numbers["Jet100U"]                  << " & " << fired_all_trigger_prescales["Jet100U"]                  << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet140U} & "            << all_trigger_numbers["Jet140U"]                  << " & " << all_trigger_prescales["Jet140U"]                  << " & " << fired_all_trigger_numbers["Jet140U"]                  << " & " << fired_all_trigger_prescales["Jet140U"]                  << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_Jet180U} & "            << all_trigger_numbers["Jet180U"]                  << " & " << all_trigger_prescales["Jet180U"]                  << " & " << fired_all_trigger_numbers["Jet180U"]                  << " & " << fired_all_trigger_prescales["Jet180U"]                  << " \\\\" << endl;   
   output_file << "\\hline" << endl;
   output_file << "Di-jet & \\texttt{HLT\\_DiJetAve15U} & "  << all_trigger_numbers["DiJetAve15U"]              << " & " << all_trigger_prescales["DiJetAve15U"]              << " & " << fired_all_trigger_numbers["DiJetAve15U"]              << " & " << fired_all_trigger_prescales["DiJetAve15U"]              << " \\\\" << endl;
   output_file << "&\\texttt{* HLT\\_DiJetAve30U} & "        << all_trigger_numbers["DiJetAve30U"]              << " & " << all_trigger_prescales["DiJetAve30U"]              << " & " << fired_all_trigger_numbers["DiJetAve30U"]              << " & " << fired_all_trigger_prescales["DiJetAve30U"]              << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_DiJetAve50U} & "        << all_trigger_numbers["DiJetAve50U"]              << " & " << all_trigger_prescales["DiJetAve50U"]              << " & " << fired_all_trigger_numbers["DiJetAve50U"]              << " & " << fired_all_trigger_prescales["DiJetAve50U"]              << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_DiJetAve70U} & "        << all_trigger_numbers["DiJetAve70U"]              << " & " << all_trigger_prescales["DiJetAve70U"]              << " & " << fired_all_trigger_numbers["DiJetAve70U"]              << " & " << fired_all_trigger_prescales["DiJetAve70U"]              << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_DiJetAve100U} & "       << all_trigger_numbers["DiJetAve100U"]             << " & " << all_trigger_prescales["DiJetAve100U"]             << " & " << fired_all_trigger_numbers["DiJetAve100U"]             << " & " << fired_all_trigger_prescales["DiJetAve100U"]             << " \\\\" << endl;   
   output_file << "&\\texttt{* HLT\\_DiJetAve140U} & "       << all_trigger_numbers["DiJetAve140U"]             << " & " << all_trigger_prescales["DiJetAve140U"]             << " & " << fired_all_trigger_numbers["DiJetAve140U"]             << " & " << fired_all_trigger_prescales["DiJetAve140U"]             << " \\\\" << endl;   
   output_file << "\\hline" << endl;
   output_file << "Quad-jet & \\texttt{HLT\\_QuadJet20U} & " << all_trigger_numbers["QuadJet20U"]               << " & " << all_trigger_prescales["QuadJet20U"]               << " & " << fired_all_trigger_numbers["QuadJet20U"]               << " & " << fired_all_trigger_prescales["QuadJet20U"]               << " \\\\" << endl;
   output_file << "& \\texttt{HLT\\_QuadJet25U} & "          << all_trigger_numbers["QuadJet25U"]               << " & " << all_trigger_prescales["QuadJet25U"]               << " & " << fired_all_trigger_numbers["QuadJet25U"]               << " & " << fired_all_trigger_prescales["QuadJet25U"]               << " \\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "$H_T$ & \\texttt{HLT\\_HT100U} & "        << all_trigger_numbers["HT100U"]                   << " & " << all_trigger_prescales["HT100U"]                   << " & " << fired_all_trigger_numbers["HT100U"]                   << " & " << fired_all_trigger_prescales["HT100U"]                   << " \\\\" << endl;
   output_file << "&\\texttt{HLT\\_HT120U} & "               << all_trigger_numbers["HT120U"]                   << " & " << all_trigger_prescales["HT120U"]                   << " & " << fired_all_trigger_numbers["HT120U"]                   << " & " << fired_all_trigger_prescales["HT120U"]                   << " \\\\" << endl;
   output_file << "&\\texttt{HLT\\_HT140U} & "               << all_trigger_numbers["HT140U"]                   << " & " << all_trigger_prescales["HT140U"]                   << " & " << fired_all_trigger_numbers["HT140U"]                   << " & " << fired_all_trigger_prescales["HT140U"]                   << " \\\\" << endl;
   output_file << "&\\texttt{HLT\\_EcalOnly\\_SumEt160} & "  << all_trigger_numbers["EcalOnly_SumEt160"]        << " & " << all_trigger_prescales["EcalOnly_SumEt160"]        << " & " << fired_all_trigger_numbers["EcalOnly_SumEt160"]        << " & " << fired_all_trigger_prescales["EcalOnly_SumEt160"]        << " \\\\" << endl;

   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;


   output_file << endl << endl << endl << endl << endl;




   output_file << "Events \t Fraction" << endl;
   output_file << "Jet Primary Dataset \t 20,022,826 \t 1.000" << endl;
   output_file << "Validated Run \t " << total_validated_events << " \t " << fixed << setprecision(3) << (total_validated_events/float(grand_total)) << endl;
   output_file << "Assigned Trigger Fired \t " << assigned_trigger_fired << " \t " << fixed << setprecision(3) << (assigned_trigger_fired/float(grand_total)) <<  endl;
   
   output_file << "Loose Jet Quality \t " << jet_quality["LOOSE"] << " \t " << fixed << setprecision(3) << (jet_quality["LOOSE"]/float(grand_total)) << endl;
   output_file << "AK5 Match \t " << ak5_match << " \t " << fixed << setprecision(3) << (ak5_match/float(grand_total)) << endl;
   output_file << "|eta| < 2.4 \t " << eta_24_cut << " \t " << fixed << setprecision(3) << (eta_24_cut/float(grand_total)) << endl;
   output_file << "Passes Soft Drop \t " << passes_softdrop[0.1] << " \t " << fixed << setprecision(3) << (passes_softdrop[0.1]/float(grand_total)) << endl;


   output_file << endl << endl << endl << endl << endl;




   
   output_file << "N_PV \t All \t All \t Filtered \t Filtered" << endl;
   
   output_file << "N_PV \t Events \t Fraction \t Events \t Fraction" << endl;
   

   output_file << "1 \t " << npv[1] << " \t " << (npv[1] / npv_total) << " \t " << filtered_npv[1] << " \t " << (filtered_npv[1] / filtered_npv_total) << endl;
   output_file << "2 \t " << npv[2] << " \t " << (npv[2] / npv_total) << " \t " << filtered_npv[2] << " \t " << (filtered_npv[2] / filtered_npv_total) << endl;
   output_file << "3 \t " << npv[3] << " \t " << (npv[3] / npv_total) << " \t " << filtered_npv[3] << " \t " << (filtered_npv[3] / filtered_npv_total) << endl;
   output_file << "4 \t " << npv[4] << " \t " << (npv[4] / npv_total) << " \t " << filtered_npv[4] << " \t " << (filtered_npv[4] / filtered_npv_total) << endl;
   output_file << "5 \t " << npv[5] << " \t " << (npv[5] / npv_total) << " \t " << filtered_npv[5] << " \t " << (filtered_npv[5] / filtered_npv_total) << endl;   
   output_file << "6 \t " << npv[6] << " \t " << (npv[6] / npv_total) << " \t " << filtered_npv[6] << " \t " << (filtered_npv[6] / filtered_npv_total) << endl;   
   output_file << "7 \t " << npv[7] << " \t " << (npv[7] / npv_total) << " \t " << filtered_npv[7] << " \t " << (filtered_npv[7] / filtered_npv_total) << endl;   
   output_file << "8 \t " << npv[8] << " \t " << (npv[8] / npv_total) << " \t " << filtered_npv[8] << " \t " << (filtered_npv[8] / filtered_npv_total) << endl;   
   output_file << "9 \t " << npv[9] << " \t " << (npv[9] / npv_total) << " \t " << filtered_npv[9] << " \t " << (filtered_npv[9] / filtered_npv_total) << endl;   
   output_file << "10 \t " << npv[10] << " \t " << (npv[10] / npv_total) << " \t " << filtered_npv[10] << " \t " << (filtered_npv[10] / filtered_npv_total) << endl;   
   output_file << "11 \t " << npv[11] << " \t " << (npv[11] / npv_total) << " \t " << filtered_npv[11] << " \t " << (filtered_npv[11] / filtered_npv_total) << endl;   
   output_file << "12 \t " << npv[12] << " \t " << (npv[12] / npv_total) << " \t " << filtered_npv[12] << " \t " << (filtered_npv[12] / filtered_npv_total) << endl;   
   output_file << "13 \t " << npv[13] << " \t " << (npv[13] / npv_total) << " \t " << filtered_npv[13] << " \t " << (filtered_npv[13] / filtered_npv_total) << endl;   
   output_file << "14 \t " << npv[14] << " \t " << (npv[14] / npv_total) << " \t " << filtered_npv[14] << " \t " << (filtered_npv[14] / filtered_npv_total) << endl;   
   output_file << ">= 15 \t " << npv_larger_than_5_total << " \t " << (npv_larger_than_5_total / npv_total) << " \t " << filtered_npv_larger_than_5_total << " \t " << (filtered_npv_larger_than_5_total / filtered_npv_total) << endl;   
   


   output_file << endl << endl << endl << endl << endl;


   output_file << "Hardest Jet p_T \t Trigger Name \t Events  \t <Prescale>" << endl;
   
   output_file << "[85, 115] GeV \t HLTJet30U \t "   << trigger_numbers["Jet30U"]     << " \t " << trigger_prescales["Jet30U"]     << endl;
   output_file << "[115, 150] GeV \t HLTJet50U \t "  << trigger_numbers["Jet50U"]     << " \t " << trigger_prescales["Jet50U"]     << endl;
   output_file << "[150, 200] GeV \t HLTJet70U  \t " << trigger_numbers["Jet70U"]     << " \t " << trigger_prescales["Jet70U"]     << endl;
   output_file << "[200, 250] GeV \t HLTJet100U \t " << trigger_numbers["Jet100U"]    << " \t " << trigger_prescales["Jet100U"]    << endl;
   
   output_file << "> 250 GeV \t Jet100U \t "         << trigger_numbers["extra_100U"] << " \t " << trigger_prescales["extra_100U"] << endl;
   output_file << "> 250 GeV \t Jet140U \t "         << trigger_numbers["Jet140U"]    << " \t " << trigger_prescales["Jet140U_"]   << endl;
   

   output_file << endl << endl << "Total number of events: " << eta_24_cut;


   output_file << endl << endl << endl;

   
	output_file << "Trigger \t Trig Present \t <Prescale> \t Trig Fired? \t <Prescale>" << endl;

	output_file << "Single-jet \t HLT_Jet15U \t "     << all_trigger_numbers["Jet15U"]                    << " \t " << all_trigger_prescales["Jet15U"]                   << " \t " << fired_all_trigger_numbers["Jet15U"]                   << " \t " << fired_all_trigger_prescales["Jet15U"]                   << endl;
	output_file << "Single-jet \t HLT_Jet15U_HNF \t " << all_trigger_numbers["Jet15U_HcalNoiseFiltered"]  << " \t " << all_trigger_prescales["Jet15U_HcalNoiseFiltered"] << " \t " << fired_all_trigger_numbers["Jet15U_HcalNoiseFiltered"] << " \t " << fired_all_trigger_prescales["Jet15U_HcalNoiseFiltered"] << endl;   
	output_file << "Single-jet \t HLT_Jet30U \t "     << all_trigger_numbers["Jet30U"]                    << " \t " << all_trigger_prescales["Jet30U"]                   << " \t " << fired_all_trigger_numbers["Jet30U"]                   << " \t " << fired_all_trigger_prescales["Jet30U"]                   << endl;   
	output_file << "Single-jet \t HLT_Jet50U \t "     << all_trigger_numbers["Jet50U"]                    << " \t " << all_trigger_prescales["Jet50U"]                   << " \t " << fired_all_trigger_numbers["Jet50U"]                   << " \t " << fired_all_trigger_prescales["Jet50U"]                   << endl;   
	output_file << "Single-jet \t HLT_Jet70U \t "     << all_trigger_numbers["Jet70U"]                    << " \t " << all_trigger_prescales["Jet70U"]                   << " \t " << fired_all_trigger_numbers["Jet70U"]                   << " \t " << fired_all_trigger_prescales["Jet70U"]                   << endl;   
	output_file << "Single-jet \t HLT_Jet100U \t "    << all_trigger_numbers["Jet100U"]                   << " \t " << all_trigger_prescales["Jet100U"]                  << " \t " << fired_all_trigger_numbers["Jet100U"]                  << " \t " << fired_all_trigger_prescales["Jet100U"]                  << endl;   
	output_file << "Single-jet \t HLT_Jet140U \t "    << all_trigger_numbers["Jet140U"]                   << " \t " << all_trigger_prescales["Jet140U"]                  << " \t " << fired_all_trigger_numbers["Jet140U"]                  << " \t " << fired_all_trigger_prescales["Jet140U"]                  << endl;   
	output_file << "Single-jet \t HLT_Jet180U \t "    << all_trigger_numbers["Jet180U"]                   << " \t " << all_trigger_prescales["Jet180U"]                  << " \t " << fired_all_trigger_numbers["Jet180U"]                  << " \t " << fired_all_trigger_prescales["Jet180U"]                  << endl;   


	output_file << "Di-jet \t HLT_DiJetAve15U \t "    << all_trigger_numbers["DiJetAve15U"]               << " \t " << all_trigger_prescales["DiJetAve15U"]              << " \t " << fired_all_trigger_numbers["DiJetAve15U"]              << " \t " << fired_all_trigger_prescales["DiJetAve15U"]              << endl;
	output_file << "Di-jet \t HLT_DiJetAve30U \t "    << all_trigger_numbers["DiJetAve30U"]               << " \t " << all_trigger_prescales["DiJetAve30U"]              << " \t " << fired_all_trigger_numbers["DiJetAve30U"]              << " \t " << fired_all_trigger_prescales["DiJetAve30U"]              << endl;   
	output_file << "Di-jet \t HLT_DiJetAve50U \t "    << all_trigger_numbers["DiJetAve50U"]               << " \t " << all_trigger_prescales["DiJetAve50U"]              << " \t " << fired_all_trigger_numbers["DiJetAve50U"]              << " \t " << fired_all_trigger_prescales["DiJetAve50U"]              << endl;   
	output_file << "Di-jet \t HLT_DiJetAve70U \t "    << all_trigger_numbers["DiJetAve70U"]               << " \t " << all_trigger_prescales["DiJetAve70U"]              << " \t " << fired_all_trigger_numbers["DiJetAve70U"]              << " \t " << fired_all_trigger_prescales["DiJetAve70U"]              << endl;   
	output_file << "Di-jet \t HLT_DiJetAve100U \t "   << all_trigger_numbers["DiJetAve100U"]              << " \t " << all_trigger_prescales["DiJetAve100U"]             << " \t " << fired_all_trigger_numbers["DiJetAve100U"]             << " \t " << fired_all_trigger_prescales["DiJetAve100U"]             << endl;   
	output_file << "Di-jet \t HLT_DiJetAve140U \t "   << all_trigger_numbers["DiJetAve140U"]              << " \t " << all_trigger_prescales["DiJetAve140U"]             << " \t " << fired_all_trigger_numbers["DiJetAve140U"]             << " \t " << fired_all_trigger_prescales["DiJetAve140U"]             << endl;   

	output_file << "Quad-jet \t HLT_QuadJet20U \t "   << all_trigger_numbers["QuadJet20U"]                << " \t " << all_trigger_prescales["QuadJet20U"]               << " \t " << fired_all_trigger_numbers["QuadJet20U"]               << " \t " << fired_all_trigger_prescales["QuadJet20U"]               << endl;
	output_file << "Quad-jet \t HLT_QuadJet25U \t "   << all_trigger_numbers["QuadJet25U"]                << " \t " << all_trigger_prescales["QuadJet25U"]               << " \t " << fired_all_trigger_numbers["QuadJet25U"]               << " \t " << fired_all_trigger_prescales["QuadJet25U"]               << endl;
	
	output_file << "H_T \t HLT_HT100U \t "            << all_trigger_numbers["HT100U"]                    << " \t " << all_trigger_prescales["HT100U"]                   << " \t " << fired_all_trigger_numbers["HT100U"]                   << " \t " << fired_all_trigger_prescales["HT100U"]                   << endl;
	output_file << "H_T \t HLT_HT120U \t "            << all_trigger_numbers["HT120U"]                    << " \t " << all_trigger_prescales["HT120U"]                   << " \t " << fired_all_trigger_numbers["HT120U"]                   << " \t " << fired_all_trigger_prescales["HT120U"]                   << endl;
	output_file << "H_T \t HLT_HT140U \t "            << all_trigger_numbers["HT140U"]                    << " \t " << all_trigger_prescales["HT140U"]                   << " \t " << fired_all_trigger_numbers["HT140U"]                   << " \t " << fired_all_trigger_prescales["HT140U"]                   << endl;
	output_file << "H_T \t HLT_EcalOnly_SumEt160 \t " << all_trigger_numbers["EcalOnly_SumEt160"]         << " \t " << all_trigger_prescales["EcalOnly_SumEt160"]        << " \t " << fired_all_trigger_numbers["EcalOnly_SumEt160"]        << " \t " << fired_all_trigger_prescales["EcalOnly_SumEt160"]        << endl;
   

   output_file << endl << endl << endl << endl << endl;






}




void get_all_files_to_process(std::vector<string> & all_files, vector<boost::filesystem::path> input_paths) {
   
   for (unsigned i = 0; i < input_paths.size(); i++) {
       boost::filesystem::path input_path = input_paths[i];
       
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
             get_all_files_to_process(all_files, { itr->path() });
          }
       }
    }
}

void count_events(MOD::Event & event_being_read, boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<int, int> & npv, boost::unordered_map<int, int> & filtered_npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut) {
   

    // second part => (?<!HLT)_([a-zA-Z][a-zA-Z][a-zA-Z\d]+)


    // First part => HLT_([A-Za-z\d]+)*


   PseudoJet trigger_jet = event_being_read.trigger_jet();
   MOD::Trigger assigned_trigger = event_being_read.assigned_trigger();
    

   // cout << assigned_trigger.name() << endl;

   // Assigned Trigger Fired. 
   if (event_being_read.assigned_trigger_fired()) {
      assigned_trigger_fired++;

      // Everybody stand back, I know regular expressions.
      boost::regex rgx_1("HLT_([A-Za-z\\d]+)*");
      boost::regex rgx_2("(?<!HLT)_([a-zA-Z][a-zA-Z][a-zA-Z\\d]+)");
      
      string assigned_trigger_name;


      boost::smatch what;
      if (boost::regex_search(assigned_trigger.name(), what, rgx_1)) {
        
        string first_part(what[1]);
        assigned_trigger_name = first_part;

        boost::smatch what_2;
        if (boost::regex_search(assigned_trigger.name(), what_2, rgx_2)) {
            string second_part(what_2[1]);
            assigned_trigger_name += second_part;
        }
      }



      if ( (assigned_trigger_name.find("100U") != std::string::npos) && (trigger_jet.pt() * trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC() > 250.0)) {
         assigned_trigger_name = "extra_100U";
      }



      auto numbers_search = trigger_numbers.find(assigned_trigger_name);


      if (numbers_search != trigger_numbers.end()) {
	      numbers_search->second++;   
	   }
	   else {
         trigger_numbers.insert(make_pair(assigned_trigger_name, 1));
	   }

      
      // string looking_for = "100U";


      // if (assigned_trigger_name.find(looking_for) != std::string::npos) {
      //   // cout << "Found." << endl;
   
      //   // cout << fixed << setprecision(8) << trigger_jet.pt() * trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC() << ", ";
      //    // cout << event_being_read.event_number() << " => " << trigger_jet.pt() * trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC() << " " << trigger_we_are_using << "....; ";
      // }


      
      

	   // Average prescale.
	   auto prescale_search = trigger_prescales.find(assigned_trigger_name);

	   if (prescale_search != trigger_prescales.end()) {

	      if (trigger_jet.pt() * trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC() > 250) {
	         if (assigned_trigger_name == "Jet140U") {
	            // Get the total number of prescales already in the hashmap.
	            int n = trigger_numbers[assigned_trigger_name];
	            float summation_x = n * trigger_prescales[assigned_trigger_name];
	            float new_mean = (summation_x + assigned_trigger.prescale()) / (n + 1);

	            trigger_prescales[assigned_trigger_name] = new_mean;      
	         }
	         else {
	            auto search_100U = trigger_numbers.find("extra_100U");
	            if (search_100U != trigger_numbers.end()) {
	               // Get the total number of prescales already in the hashmap.
	               int n = trigger_numbers[assigned_trigger_name];
	               float summation_x = n * trigger_prescales[assigned_trigger_name];
	               float new_mean = (summation_x + assigned_trigger.prescale()) / (n + 1);

	               trigger_prescales["extra_100U"] = new_mean;      
	            }
	            else {
	               trigger_prescales.insert(make_pair("extra_100U", assigned_trigger.prescale()));
	            }
	               
	         }
	      }
	      else {
	         // Get the total number of prescales already in the hashmap.
	         int n = trigger_numbers[assigned_trigger_name];
	         float summation_x = n * trigger_prescales[assigned_trigger_name];
	         float new_mean = (summation_x + assigned_trigger.prescale()) / (n + 1);

	         trigger_prescales[assigned_trigger_name] = new_mean;      
	      }
	      
	   }
	   else {
	      if (trigger_jet.pt() * trigger_jet.user_info<MOD::InfoCalibratedJet>().JEC() > 250) {
	         if (assigned_trigger_name == "Jet140U") {
	            trigger_prescales.insert(make_pair(assigned_trigger_name, assigned_trigger.prescale()));
	         }
	         else {
	            trigger_prescales.insert(make_pair("extra_100U", assigned_trigger.prescale()));
	         }     
	      }
	      else {
	         trigger_prescales.insert(make_pair(assigned_trigger_name, assigned_trigger.prescale()));
	      }
	   }




      // Loose jet quality.
      
      // First, find the trigger jet.
      
      string quality_string;

      if (trigger_jet.has_user_info()) {
         int quality = trigger_jet.user_info<MOD::InfoCalibratedJet>().jet_quality();

         // cout << quality << ", ";

         if (quality == 0) {
            quality_string = "FAILED";
         }
         else if (quality > 0) {
         	quality_string = "LOOSE";
         }
         else {
         	quality_string = "ERROR";
         	// cout << "ERROR" << quality << endl;
         }
            
         
         auto search = jet_quality.find(quality_string);

         if (search != jet_quality.end()) {
            search->second++;            
         }
         else {
            jet_quality.insert(make_pair(quality_string, 1));
         }
      }


      if (quality_string == "LOOSE") {
         // AK5 match.
         if (event_being_read.is_trigger_jet_matched()) {
            ak5_match++;
            
            // Eta cut.
            if (abs(event_being_read.hardest_jet().eta()) < 2.4) {
               eta_24_cut++;





               // NPV.
               // NPV
               auto search = filtered_npv.find(event_being_read.condition().npv());

               if (search != filtered_npv.end())
                  search->second++;
               else
                  filtered_npv.insert(make_pair(event_being_read.condition().npv(), 1));

               // // SoftDrop

               // // Passes SoftDrop
   
               JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);

               if (event_being_read.hardest_jet().has_structure()) {

                  Selector pT_1_GeV_selector = SelectorPtMin(1.0);
                  vector<PseudoJet> hardest_jet_constituents = pT_1_GeV_selector(event_being_read.hardest_jet().constituents());

                  ClusterSequence cs_uncorrected_jet_no_softkiller = ClusterSequence(hardest_jet_constituents, jet_def_cambridge);   

                  PseudoJet uncorrected_hardest_jet_no_softkiller = cs_uncorrected_jet_no_softkiller.inclusive_jets()[0];
                  
                  if (uncorrected_hardest_jet_no_softkiller.has_structure()) {
                     SoftDrop soft_drop_005(0.0, 0.05);
                     SoftDrop soft_drop_01(0.0, 0.1);
                     SoftDrop soft_drop_02(0.0, 0.2);

                     PseudoJet soft_drop_005_jet = soft_drop_005(uncorrected_hardest_jet_no_softkiller);
                     PseudoJet soft_drop_01_jet = soft_drop_01(uncorrected_hardest_jet_no_softkiller);
                     PseudoJet soft_drop_02_jet = soft_drop_02(uncorrected_hardest_jet_no_softkiller);

                     double zg_005 = soft_drop_005_jet.structure_of<SoftDrop>().symmetry();
                     double zg_01 = soft_drop_01_jet.structure_of<SoftDrop>().symmetry();
                     double zg_02 = soft_drop_02_jet.structure_of<SoftDrop>().symmetry();

                     if (zg_005 > 0.05) {
                        auto search = passes_softdrop.find(0.05);

                        if (search != passes_softdrop.end())
                           search->second++;
                        else
                           passes_softdrop.insert(make_pair(0.05, 1));
                     }

                     if (zg_01 > 0.1) {
                        auto search = passes_softdrop.find(0.1);

                        if (search != passes_softdrop.end())
                           search->second++;
                        else
                           passes_softdrop.insert(make_pair(0.1, 1));
                     }

                     if (zg_02 > 0.2) {
                        auto search = passes_softdrop.find(0.2);

                        if (search != passes_softdrop.end())
                           search->second++;
                        else
                           passes_softdrop.insert(make_pair(0.2, 1));
                     }
                  }
                  
               }




            }
         }
         else {
         	cout << "Found an event in loose for which there is no AK5 match." << endl;
         }


      }


   }
      
   // NPV
   auto search = npv.find(event_being_read.condition().npv());

   if (search != npv.end())
      search->second++;
   else
      npv.insert(make_pair(event_being_read.condition().npv(), 1));  
   
	

}










void count_triggers(MOD::Event & event_being_read, boost::unordered_map<string, int> & all_trigger_numbers, boost::unordered_map<string, float> & all_trigger_prescales, boost::unordered_map<string, int> & all_fired_trigger_numbers, boost::unordered_map<string, float> & all_fired_trigger_prescales) {
   
   vector<MOD::Trigger> all_triggers = event_being_read.triggers();

   for (unsigned i=0; i < all_triggers.size(); i++) {
      MOD::Trigger trigger = all_triggers[i];

      string trigger_name;

      // Everybody stand back, I know regular expressions.
      boost::regex rgx_1("HLT_([A-Za-z0-9]*)");
      boost::regex rgx_2("(?<!HLT)_([a-zA-Z][a-zA-Z][a-zA-Z0-9]+)");
      
    
      string given_trigger_name = trigger.name();

      string first_part = "";
      string second_part = "";
      
      boost::smatch what;
      if (boost::regex_search(given_trigger_name, what, rgx_1)) {
        
        first_part = what[1];
        
        boost::smatch what_2;
        if (boost::regex_search(given_trigger_name, what_2, rgx_2)) {
            second_part = "_" + what_2[1];
        }
      }


      trigger_name = first_part + second_part;
      


      // First present triggers.
      auto numbers_search = all_trigger_numbers.find(trigger_name);

      if (numbers_search != all_trigger_numbers.end()) {
         numbers_search->second++;
      }
      else {
         all_trigger_numbers.insert(make_pair(trigger_name, 1));
      }

      // Average prescale.
      auto prescale_search = all_trigger_prescales.find(trigger_name);

      if (prescale_search != all_trigger_prescales.end()) {
         // Get the total number of prescales already in the hashmap.
         int n = all_trigger_numbers[trigger_name];
         float summation_x = n * all_trigger_prescales[trigger_name];
         float new_mean = (summation_x + trigger.prescale()) / (n + 1);

         all_trigger_prescales[trigger_name] = new_mean;   
      }
      else {
         all_trigger_prescales.insert(make_pair(trigger_name, trigger.prescale()));
      }


      // Next, fired triggers.
      if (trigger.fired()) {
         auto numbers_search = all_fired_trigger_numbers.find(trigger_name);

         if (numbers_search != all_fired_trigger_numbers.end()) {
            numbers_search->second++;
         }
         else {
            all_fired_trigger_numbers.insert(make_pair(trigger_name, 1));
         }

         // Average prescale.
         auto prescale_search = all_fired_trigger_prescales.find(trigger_name);

         if (prescale_search != all_fired_trigger_prescales.end()) {
            // Get the total number of prescales already in the hashmap.
            int n = all_fired_trigger_numbers[trigger_name];
            float summation_x = n * all_fired_trigger_prescales[trigger_name];
            float new_mean = (summation_x + trigger.prescale()) / (n + 1);

            all_fired_trigger_prescales[trigger_name] = new_mean;   
         }
         else {
            all_fired_trigger_prescales.insert(make_pair(trigger_name, trigger.prescale()));
         }
      }
      
   }

   
}
