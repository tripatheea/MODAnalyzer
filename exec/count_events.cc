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


void count_events(MOD::Event & event_being_read, boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<int, int> & npv, boost::unordered_map<int, int> & filtered_npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut);
void get_all_files_to_process(std::vector<string> & all_files, boost::filesystem::path input_path);

void write_stats(boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<int, int> & npv, boost::unordered_map<int, int> & filtered_npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut, int & total_validated_events, int & grand_total);

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
   
   cout << endl << endl << "Starting count_events with the following given arguments: " << endl << endl;
   cout << "Input Path   : " << argv[1] << endl;
   cout << "Use LaTeX    : " << use_latex << endl;
   cout << "No. of Files : " << number_of_files_to_process << endl;
   cout << "No. of Events: " << number_of_events_to_process << endl;
   
   ofstream output_file("stats.dat", ios::out);

   // Recursively collect all filenames to process.

   vector<string> all_filenames;
   get_all_files_to_process(all_filenames, input_path);

   // Setup data structures to hold all counts. 
   boost::unordered_map<int, int> npv;
   boost::unordered_map<int, int> filtered_npv;
   boost::unordered_map<string, int> jet_quality;
   boost::unordered_map<float, int> passes_softdrop;

   boost::unordered_map<string, int> trigger_numbers;
   boost::unordered_map<string, float> trigger_prescales;

   // int total_files = all_filenames.size();
   int total_files = all_filenames.size();
   int total_validated_events = 0;
   int assigned_trigger_fired = 0;
   int ak5_match = 0;
   int eta_24_cut = 0;
   

   int grand_total = 20022826;
   // int grand_total = 20022;

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
            
            count_events(event_being_read, trigger_numbers, trigger_prescales, npv, filtered_npv, jet_quality, passes_softdrop, assigned_trigger_fired, ak5_match, eta_24_cut);

            event_being_read = MOD::Event();
            event_serial_number++;
            total_validated_events++;
         }   
      }
      else {
         break;
      }

      write_stats(trigger_numbers, trigger_prescales, npv, filtered_npv, jet_quality, passes_softdrop, assigned_trigger_fired, ak5_match, eta_24_cut, total_validated_events, grand_total);
      
   }

   
   cout << endl << "Everything done. Printing summary." << endl << endl;

   cout << endl << "==================================================================" << endl << endl;



   cout << "Assigned Trigger Fired: " << assigned_trigger_fired << endl << endl;
   cout << "AK5 Match: " << ak5_match << endl << endl;
   cout << "|eta| < 2.4: " << eta_24_cut << endl << endl;

   cout << "====== NPV ======" << endl;
   for (auto kv : npv) {
      cout << kv.first << " => " << kv.second << endl;
   }
   cout << endl << endl;

   cout << "====== Jet Quality ======" << endl;
   for (auto kv : jet_quality) {
      cout << kv.first << " => " << kv.second << endl;
   }
   cout << endl << endl;

   cout << "====== Total ======" << endl;
   cout << "Total Validated Events: " << total_validated_events << endl;   


   if (use_latex) {


      cout << endl << endl << "================ LaTeX =============" << endl << endl;

   
      cout << "\\hline" << endl;
      cout << "\\hline" << endl;
      cout << "& Events & Fraction\\\\" << endl;
      cout << "\\hline" << endl;
      cout << "Jet Primary Dataset & 20,022,826 & 1.000 \\\\" << endl;
      cout << "Validated Run & " << total_validated_events << " & " << fixed << setprecision(3) << (total_validated_events/float(grand_total)) << " \\\\" << endl;
      cout << "Loose Jet Quality (\\Tab{tab:jet_quality}) & " << jet_quality["LOOSE"] << " & " << fixed << setprecision(3) << (jet_quality["LOOSE"]/float(grand_total)) << " \\\\" << endl;
      cout << "Assigned Trigger Fired (\\Tab{tab:trigger_table}) & " << assigned_trigger_fired << " & " << fixed << setprecision(3) << (assigned_trigger_fired/float(grand_total)) << " \\\\" << endl;
      cout << "\\hline" << endl;
      cout << "AK5 Match & " << ak5_match << " & " << fixed << setprecision(3) << (ak5_match/float(grand_total)) << "\\\\" << endl;
      cout << "$\\left| \\eta \\right| < 2.4$ & " << eta_24_cut << " & " << fixed << setprecision(3) << (eta_24_cut/float(grand_total)) << " \\\\" << endl;
      cout << "\\hline" << endl;
      cout << "Passes Soft Drop ($z_g > z_{\\rm cut}$) & " << passes_softdrop[0.1] << " & " << fixed << setprecision(3) << (passes_softdrop[0.1]/float(grand_total)) << " \\\\" << endl;
      cout << "\\hline" << endl;
      cout << "\\hline" << endl;

      
      write_stats(trigger_numbers, trigger_prescales, npv, filtered_npv, jet_quality, passes_softdrop, assigned_trigger_fired, ak5_match, eta_24_cut, total_validated_events, grand_total);
         

      
   }

   cout << endl << endl;
   

}


void write_stats(boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<int, int> & npv, boost::unordered_map<int, int> & filtered_npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut, int & total_validated_events, int & grand_total) {
   // cout << endl << endl << "================ LaTeX =============" << endl << endl;

   ofstream output_file("./stats.dat", ios::out);

   // ostringstream LaTeX_stream;

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
      if (kv.first >= 6)
         npv_larger_than_5_total += kv.second;
   }

   float filtered_npv_total = 0;
   float filtered_npv_larger_than_5_total = 0;
   for (auto kv : filtered_npv) {
      // output_file << kv.first << " => " << kv.second << endl;
      filtered_npv_total += kv.second;
      if (kv.first >= 6)
         filtered_npv_larger_than_5_total += kv.second;
   }


   // output_file << endl << endl;

   

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
   output_file << "$\\geq 6$ & " << npv_larger_than_5_total << " & " << (npv_larger_than_5_total / npv_total) << " & " << filtered_npv_larger_than_5_total << " & " << (filtered_npv_larger_than_5_total / filtered_npv_total) << " \\\\" << endl;
   
   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;

   // cout << LaTeX_stream.str() << endl;   

   output_file << endl << endl << endl << endl << endl;

   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;
   output_file << "Hardest Jet $p_T$ & Trigger Name & Events  & $\\langle$Prescale$\\rangle$ \\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "$[85, 115]~\\GeV$ & \\texttt{HLT\\_Jet30U} & " << trigger_numbers["HLT_Jet30U"] + trigger_numbers["HLT_Jet30U_v3"] << " & " << (trigger_prescales["HLT_Jet30U"] * trigger_numbers["HLT_Jet30U"] +  trigger_prescales["HLT_Jet30U_v3"] * trigger_numbers["HLT_Jet30U_v3"]) / (trigger_numbers["HLT_Jet30U"] + trigger_numbers["HLT_Jet30U_v3"]) << " \\\\" << endl;
   output_file << "$[115, 150]~\\GeV$ & \\texttt{HLT\\_Jet50U} & " << trigger_numbers["HLT_Jet50U"] + trigger_numbers["HLT_Jet50U_v3"] << " & " << (trigger_prescales["HLT_Jet50U"] * trigger_numbers["HLT_Jet50U"] +  trigger_prescales["HLT_Jet50U_v3"] * trigger_numbers["HLT_Jet50U_v3"]) / (trigger_numbers["HLT_Jet50U"] + trigger_numbers["HLT_Jet50U_v3"]) << " \\\\" << endl;
   output_file << "$[150, 200]~\\GeV$ & \\texttt{HLT\\_Jet70U}  & " << trigger_numbers["HLT_Jet70U"] + trigger_numbers["HLT_Jet70U_v2"] + trigger_numbers["HLT_Jet70U_v3"] << " & " << (trigger_prescales["HLT_Jet70U"] * trigger_numbers["HLT_Jet70U"] + trigger_prescales["HLT_Jet70U_v2"] * trigger_numbers["HLT_Jet70U_v2"] +  trigger_prescales["HLT_Jet70U_v3"] * trigger_numbers["HLT_Jet70U_v3"]) / (trigger_numbers["HLT_Jet70U"] + trigger_numbers["HLT_Jet70U_v2"] + trigger_numbers["HLT_Jet70U_v3"]) << " \\\\" << endl;
   output_file << "$[200, 250]~\\GeV$ & \\texttt{HLT\\_Jet100U} & " << trigger_numbers["HLT_Jet100U"] + trigger_numbers["HLT_Jet100U_v2"] + trigger_numbers["HLT_Jet100U_v3"] << " & " << (trigger_prescales["HLT_Jet100U"] * trigger_numbers["HLT_Jet100U"] + trigger_prescales["HLT_Jet100U_v2"] * trigger_numbers["HLT_Jet100U_v2"] +  trigger_prescales["HLT_Jet100U_v3"] * trigger_numbers["HLT_Jet100U_v3"]) / (trigger_numbers["HLT_Jet100U"] + trigger_numbers["HLT_Jet100U_v2"] + trigger_numbers["HLT_Jet100U_v3"]) <<  "\\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "\\multirow{2}{*}{$> 250~\\GeV$} & \\texttt{HLT\\_Jet100U}   & " << trigger_numbers["extra_100U"] << " & " << trigger_prescales["extra_100U"] << "\\\\" << endl;
   output_file << "& \\texttt{HLT\\_Jet140U} & " << trigger_numbers["HLT_Jet140U_v1"] + trigger_numbers["HLT_Jet140U_v3"] << " & " << (trigger_prescales["HLT_Jet140U_v1"] * trigger_numbers["HLT_Jet140U_v1"] +  trigger_prescales["HLT_Jet140U_v3"] * trigger_numbers["HLT_Jet140U_v3"]) / (trigger_numbers["HLT_Jet140U_v1"] + trigger_numbers["HLT_Jet140U_v3"]) << "\\\\" << endl;
   output_file << "\\hline" << endl;
   output_file << "\\hline" << endl;

   // output_file << LaTeX_stream;      
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

void count_events(MOD::Event & event_being_read, boost::unordered_map<string, int> & trigger_numbers, boost::unordered_map<string, float> & trigger_prescales, boost::unordered_map<int, int> & npv, boost::unordered_map<int, int> & filtered_npv, boost::unordered_map<string, int> & jet_quality, boost::unordered_map<float, int> & passes_softdrop, int & assigned_trigger_fired, int & ak5_match, int & eta_24_cut) {
   
   PseudoJet trigger_jet = event_being_read.trigger_jet();
   MOD::Trigger assigned_trigger = event_being_read.assigned_trigger();
   auto numbers_search = trigger_numbers.find(assigned_trigger.name());

   // if ((trigger_jet.pt() > 250)) {
      // cout << assigned_trigger.name() << endl;
   // }

   if (numbers_search != trigger_numbers.end()) {

      if (trigger_jet.pt() > 250) {
         if ((assigned_trigger.name() == "HLT_Jet140U_v1") || (assigned_trigger.name() == "HLT_Jet140U_v3")) {
            numbers_search->second++;      
         }
         else {
            auto search_100U = trigger_numbers.find("extra_100U");
            if (search_100U != trigger_numbers.end())
               search_100U->second++;
            else
               trigger_numbers.insert(make_pair("extra_100U", 1));
         }
      }
      else {
         numbers_search->second++;   
      }
      
   }
   else {
      if (trigger_jet.pt() > 250) {
         if ((assigned_trigger.name() == "HLT_Jet140U_v1") || (assigned_trigger.name() == "HLT_Jet140U_v3")) {
            trigger_numbers.insert(make_pair(assigned_trigger.name(), 1));
         }
         else {
            trigger_numbers.insert(make_pair("extra_100U", 1));
         }     
      }
      else {
         trigger_numbers.insert(make_pair(assigned_trigger.name(), 1));
      }
   }

   // Average prescale.
   auto prescale_search = trigger_prescales.find(assigned_trigger.name());

   if (prescale_search != trigger_prescales.end()) {

      if (trigger_jet.pt() > 250) {
         if ((assigned_trigger.name() == "HLT_Jet140U_v1") || (assigned_trigger.name() == "HLT_Jet140U_v3")) {
            // Get the total number of prescales already in the hashmap.
            int n = trigger_numbers[assigned_trigger.name()];
            float summation_x = n * trigger_prescales[assigned_trigger.name()];
            float new_mean = (summation_x + assigned_trigger.prescale()) / (n + 1);

            trigger_prescales[assigned_trigger.name()] = new_mean;      
         }
         else {
            auto search_100U = trigger_numbers.find("extra_100U");
            if (search_100U != trigger_numbers.end()) {
               // Get the total number of prescales already in the hashmap.
               int n = trigger_numbers[assigned_trigger.name()];
               float summation_x = n * trigger_prescales[assigned_trigger.name()];
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
         int n = trigger_numbers[assigned_trigger.name()];
         float summation_x = n * trigger_prescales[assigned_trigger.name()];
         float new_mean = (summation_x + assigned_trigger.prescale()) / (n + 1);

         trigger_prescales[assigned_trigger.name()] = new_mean;      
      }
      
   }
   else {
      if (trigger_jet.pt() > 250) {
         if ((assigned_trigger.name() == "HLT_Jet140U_v1") || (assigned_trigger.name() == "HLT_Jet140U_v3")) {
            trigger_prescales.insert(make_pair(assigned_trigger.name(), assigned_trigger.prescale()));
         }
         else {
            trigger_prescales.insert(make_pair("extra_100U", assigned_trigger.prescale()));
         }     
      }
      else {
         trigger_prescales.insert(make_pair(assigned_trigger.name(), assigned_trigger.prescale()));
      }
   }



   // Assigned Trigger Fired. 
   if (event_being_read.assigned_trigger_fired()) {
      assigned_trigger_fired++;

      // Loose jet quality.
      
      // First, find the trigger jet.
      PseudoJet trigger_jet = event_being_read.trigger_jet();

      string quality_string;

      if (trigger_jet.has_user_info()) {
         int quality = trigger_jet.user_info<MOD::InfoCalibratedJet>().jet_quality();

         if (quality == 0)
            quality_string = "FAILED";
         else if (quality == 1)
            quality_string = "LOOSE";
         else if (quality == 2)
            quality_string = "MEDIUM";
         else if (quality == 3)
            quality_string = "TIGHT";
         else
            quality_string = "ERROR";
         
         auto search = jet_quality.find(quality_string);

         if (search != jet_quality.end()) {
            
            search->second++;

            // If an event passes "Tight", it also passes "Medium" and "Loose". Implement that.

            if (quality_string == "TIGHT") {
               auto medium_search = jet_quality.find("MEDIUM"); if (medium_search != jet_quality.end()) { jet_quality.find("MEDIUM")->second++; } else { jet_quality.insert(make_pair("MEDIUM", 1)); }
               auto loose_search = jet_quality.find("LOOSE"); if (loose_search != jet_quality.end()) { jet_quality.find("LOOSE")->second++; } else { jet_quality.insert(make_pair("LOOSE", 1)); }

               if (medium_search != jet_quality.end())
                  medium_search++;
               else
                  jet_quality.insert(make_pair("MEDIUM", 1));

               if (loose_search != jet_quality.end())
                  loose_search++;
               else
                  jet_quality.insert(make_pair("LOOSE", 1));
            }
            else if (quality_string == "MEDIUM") {
               auto loose_search = jet_quality.find("LOOSE"); if (loose_search != jet_quality.end()) { jet_quality.find("LOOSE")->second++; } else { jet_quality.insert(make_pair("LOOSE", 1)); }

               if (loose_search != jet_quality.end())
                  loose_search++;
               else
                  jet_quality.insert(make_pair("LOOSE", 1));
            }

         }
         else {
            jet_quality.insert(make_pair(quality_string, 1));
         }
      }

      if ((quality_string == "LOOSE") || (quality_string == "MEDIUM") || (quality_string == "TIGHT")) {
         // AK5 match.
         if (event_being_read.trigger_jet_is_matched()) {
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

               // SoftDrop

               // Passes SoftDrop
   
               JetDefinition jet_def_cambridge(cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);

               if (event_being_read.hardest_jet().has_structure()) {
                  ClusterSequence cs_uncorrected_jet_no_softkiller = ClusterSequence(event_being_read.hardest_jet().constituents(), jet_def_cambridge);   

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
      }

   }
      
   // NPV
   auto search = npv.find(event_being_read.condition().npv());

   if (search != npv.end())
      search->second++;
   else
      npv.insert(make_pair(event_being_read.condition().npv(), 1));  
   
}
