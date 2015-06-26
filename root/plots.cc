#include <iostream>
#include <fstream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TVectorT.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TLegend.h"

#include "TMath.h"

#include "THStack.h"

#include <unordered_map>

#include <sstream>
#include <vector>


using namespace std;


void n_tilde_against_jet_multiplicity();
void hardest_pt_corresponding_triggers();
void fix_cone_radius_sweep_pt_cut();
void fix_pt_cut_sweep_cone_radius();
void corrected_ak5_spectrum();
void invariant_mass();
void zg_plots();


void plots() {
  hardest_pt_corresponding_triggers();
  // corrected_ak5_spectrum();
  // zg_plots();
}


void n_tilde_against_jet_multiplicity() {
  ifstream infile("../data/CMS_JetSample_analyzed.dat");

  TFile * rootFile_;
  TTree * multiplicityTree_;

  string tag, trigger_name;
  int event_number, run_number, antikt_jets_size, prescale;
  double n_tilde, cone_radius, pt_cut, hardest_pt_ak5, hardest_pt_ak7;
  bool fired;

  THStack *hs = new THStack("Fractional Jet Multiplicity", "Fractional Jet Multiplicity (pt_cut = 50.0 GeV, R = 0.5)");

  vector<TH1F * > N_tildes = vector<TH1F *>();
  EColor colors[7] = {kRed, kBlue, kGreen, kYellow, kMagenta, kOrange, kCyan};
  const char * antikT_labels[7] = {"antikT Jet Size = 0", "antikT Jet Size = 1", "antikT Jet Size = 2", "antikT Jet Size = 3", "antikT Jet Size = 4", "antikT Jet Size = 5", "antikT Jet Size = 6"};
  
  gStyle->SetOptStat(false);

  for (int i = 0; i < 7; i++) {
    TH1F * N_tilde_temp = new TH1F("", "", 50, -0.5, 6.0);
    N_tildes.push_back(N_tilde_temp);
  }

  string line;
  while(getline(infile, line)) {
    istringstream iss(line);
    
    iss >> tag >> event_number >> run_number >> n_tilde >> antikt_jets_size >> trigger_name >> fired >> prescale >> cone_radius >> pt_cut >> hardest_pt_ak5 >> hardest_pt_ak7;

    if (tag != "#") {
      if ((fired) && (pt_cut == 50) && (cone_radius = 0.50)) {
        N_tildes[antikt_jets_size]->Fill(n_tilde, prescale);
      }
    }
  }    


  TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

  for(int i = 0; i < 6; i++) {
    N_tildes[i]->SetFillColorAlpha(colors[i], 0.5);
    N_tildes[i]->SetMarkerStyle(21);
    N_tildes[i]->SetMarkerColor(colors[i]);
    hs->Add(N_tildes[i]);
    legend->AddEntry(N_tildes[i], antikT_labels[i]);
  }
  
  TCanvas *c2e = new TCanvas("c2e", "c2e", 600, 400);

  c2e->BuildLegend();

  gPad->SetLogy();
  
  hs->SetMaximum(10e6);
  hs->SetMinimum(10e-3);

  hs->Draw();
  hs->GetHistogram()->GetXaxis()->SetTitle("Fractional Jet Multiplicity (Jets Without Jets)");
  hs->GetHistogram()->GetXaxis()->CenterTitle();
  
  legend->Draw();
}



void hardest_pt_corresponding_triggers() {
  ifstream infile("/media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat");


  TFile * rootFile_;
  TTree * multiplicityTree_;

  string tag, trigger_name;
  int event_number, run_number, prescale;
  double uncorrected_pt, corrected_pt, zg_05, zg_1, zg_2;

  
  unordered_map<double, TH1F * > hardest_pts;

  EColor colors[5] = {kRed, kBlue, kGreen, kYellow, kMagenta};
  const char * trigger_labels[5] = {"HLT_Jet70U", "HLT_Jet50U", "HLT_Jet30U", "HLT_Jet15U", "HLT_L1Jet6U"};
  
  gStyle->SetOptStat(false);

  for (unsigned int i = 0; i < 5; i++) {
    TH1F * pt_temp = new TH1F("a", "", 50, -50, 1000);
    hardest_pts[i] = pt_temp;

  }
  
  string line;
  while(getline(infile, line)) {
    istringstream iss(line);
    
    iss >> tag >> event_number >> run_number >> uncorrected_pt >> corrected_pt >> prescale >> trigger_name >> zg_05 >> zg_1 >> zg_2;

    if (tag != "#") {
        int trigger_index = std::distance(trigger_labels, std::find(trigger_labels, trigger_labels + 6, trigger_name));
        hardest_pts[trigger_index]->Fill(corrected_pt, prescale);
    }
  }


  TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

  for(unsigned int i = 0; i < 5; i++) {
    hardest_pts[i]->SetFillColorAlpha(colors[i], 0.5);
    hardest_pts[i]->SetMarkerStyle(21);
    hardest_pts[i]->SetMarkerColor(colors[i]);
    
    legend->AddEntry(hardest_pts[i], trigger_labels[i]);

    if (i == 0) {
        hardest_pts[i]->Draw("E");
        hardest_pts[i]->GetXaxis()->SetTitle("Hardest Jet pT");
        hardest_pts[i]->GetXaxis()->CenterTitle();
        hardest_pts[i]->GetYaxis()->SetRangeUser(10e-2, 10e6);
    }
    else {
      hardest_pts[i]->Draw("same E");            
    }
  }

  gPad->SetLogy();
  
  legend->Draw();

  gPad->Print("hardest_pt_corresponding_triggers.pdf");
}




void fix_cone_radius_sweep_pt_cut() {

  // Fix R, sweep across pt_cut.

  double fixed_cone_radius = 0.5;

  ifstream infile("../data/CMS_JetSample_analyzed.dat");

  TFile * rootFile_;
  TTree * multiplicityTree_;
  
  string tag, trigger_name;
  int event_number, run_number, antikt_jets_size, prescale;
  double n_tilde, cone_radius, pt_cut, hardest_pt_ak5, hardest_pt_ak7;
  bool fired;

  THStack *hs = new THStack("Fractional Jet Multiplicity", "Fractional Jet Multiplicity (R = 0.5)");

  unordered_map<double, TH1F * > pt_cuts_map;
  vector<int> pt_cuts = {50, 80, 110};

  EColor colors[3] = {kRed, kBlue, kGreen};
  const char * pt_cuts_labels[3] = {"pT_cut = 50 GeV", "pT_cut = 80 GeV", "pT_cut = 110 GeV"};
  
  gStyle->SetOptStat(false);

  for (unsigned int i = 0; i < pt_cuts.size(); i++) {
    TH1F * pt_cut_temp = new TH1F("", "", 50, -0.5, 6.0);
    pt_cuts_map[pt_cuts[i]] = pt_cut_temp;
  }


  string line;
  while(getline(infile, line)) {
    istringstream iss(line);
    
    iss >> tag >> event_number >> run_number >> n_tilde >> antikt_jets_size >> trigger_name >> fired >> prescale >> cone_radius >> pt_cut >> hardest_pt_ak5 >> hardest_pt_ak7;

    if (tag != "#") {
      if ((fired) && (cone_radius == fixed_cone_radius)) {
        pt_cuts_map[pt_cut]->Fill(n_tilde, prescale);
      }
    }
  }


  TCanvas *cst = new TCanvas("cst","N_tilde and AntikT", 1000, 600);

  TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

  for(unsigned int i = 0; i < pt_cuts.size(); i++) {
    pt_cuts_map[pt_cuts[i]]->SetFillColorAlpha(colors[i], 0.5);
    pt_cuts_map[pt_cuts[i]]->SetMarkerStyle(21);
    pt_cuts_map[pt_cuts[i]]->SetMarkerColor(colors[i]);
    
    hs->Add(pt_cuts_map[pt_cuts[i]]);
  
    legend->AddEntry(pt_cuts_map[pt_cuts[i]], pt_cuts_labels[i]);
  }
  
  cst->BuildLegend();

  gPad->SetLogy();
  
  hs->SetMaximum(10e6);
  hs->SetMinimum(10e-3);

  hs->Draw();

  hs->GetHistogram()->GetXaxis()->SetTitle("N_tilde");
  hs->GetHistogram()->GetXaxis()->CenterTitle();

  legend->Draw();

}





void fix_pt_cut_sweep_cone_radius() {

  // Fix pt_cut, sweep across R.
  int fixed_pt_cut = 80;

  ifstream infile("../data/CMS_JetSample_analyzed.dat");

  TFile * rootFile_;
  TTree * multiplicityTree_;
  
  string tag, trigger_name;
  int event_number, run_number, antikt_jets_size, prescale;
  double n_tilde, cone_radius, pt_cut, hardest_pt_ak5, hardest_pt_ak7;
  bool fired;

  THStack *hs = new THStack("Fractional Jet Multiplicity", "Fractional Jet Multiplicity (pt_cut = 80 GeV)");

  unordered_map<double, TH1F * > cone_radii_map;
  vector<double> cone_radii = {0.3, 0.5, 0.7};

  EColor colors[3] = {kRed, kBlue, kGreen};
  const char * cone_radii_labels[3] = {"R = 0.3", "R = 0.5", "R = 0.7"};
  
  gStyle->SetOptStat(false);

  for (unsigned int i = 0; i < cone_radii.size(); i++) {
    TH1F * cone_radius_temp = new TH1F("", "", 50, -0.5, 6.0);
    cone_radii_map[cone_radii[i]] = cone_radius_temp;
  }


  string line;
  while(getline(infile, line)) {
    istringstream iss(line);
    
    iss >> tag >> event_number >> run_number >> n_tilde >> antikt_jets_size >> trigger_name >> fired >> prescale >> cone_radius >> pt_cut >> hardest_pt_ak5 >> hardest_pt_ak7;
    
    if (tag != "#") {
      if ((fired) && (pt_cut == fixed_pt_cut)) {
        cone_radii_map[cone_radius]->Fill(n_tilde, prescale);
      }
    }
  }


  TCanvas *cst = new TCanvas("cst","N_tilde and AntikT", 1000, 600);

  TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

  for(unsigned int i = 0; i < cone_radii.size(); i++) {
    cone_radii_map[cone_radii[i]]->SetFillColorAlpha(colors[i], 0.5);
    cone_radii_map[cone_radii[i]]->SetMarkerStyle(21);
    cone_radii_map[cone_radii[i]]->SetMarkerColor(colors[i]);
    
    hs->Add(cone_radii_map[cone_radii[i]]);
  
    legend->AddEntry(cone_radii_map[cone_radii[i]], cone_radii_labels[i]);
  }

  cst->BuildLegend();
  
  gPad->SetLogy();
  
  hs->SetMaximum(10e6);
  hs->SetMinimum(10e-3);

  hs->Draw();
  
  hs->GetHistogram()->GetXaxis()->SetTitle("N_tilde");
  hs->GetHistogram()->GetXaxis()->CenterTitle();

  legend->Draw();

}



void corrected_ak5_spectrum() {


  double pt_cut = 50.00;


  ifstream infile("/media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat");

  TFile * rootFile_;
  TTree * multiplicityTree_;

  
  string tag, trigger_name;
  double uncorrected_pt, corrected_pt, zg_05, zg_1, zg_2;
  int prescale, event_number, run_number;

  bool fired;

  THStack *hs = new THStack("Corrected AK5 Spectrum", "Corrected AK5 Spectrum");


  unordered_map<std::string, TH1F * > cone_radii_map;
  
  EColor colors[2] = {kRed, kGreen};
  
  gStyle->SetOptStat(false);



  cone_radii_map["Corrected"] = new TH1F("", "", 50, -20.0, 1000.0);
  cone_radii_map["Uncorrected"] = new TH1F("", "", 50, -20.0, 1000.0);



  string line;
  while(getline(infile, line)) {
    istringstream iss(line);

    iss >> tag >> event_number >> run_number >> uncorrected_pt >> corrected_pt >> prescale >> trigger_name >> zg_05 >> zg_1 >> zg_2;
    
    if (tag != "#") {
      if (uncorrected_pt > pt_cut)
        cone_radii_map["Uncorrected"]->Fill(uncorrected_pt, prescale);
      if (corrected_pt > pt_cut)
        cone_radii_map["Corrected"]->Fill(corrected_pt, prescale);

    }
  }

  TCanvas *cst = new TCanvas("cst","Corrected AK5 Spectrum", 1000, 600);

  TLegend * legend = new TLegend(0.7, 0.9, 0.9, 0.8);


  cone_radii_map["Corrected"]->SetFillColorAlpha(kRed, 0.5);
  cone_radii_map["Corrected"]->SetMarkerStyle(21);
  cone_radii_map["Corrected"]->SetMarkerColor(kRed);



  cone_radii_map["Uncorrected"]->SetFillColorAlpha(kGreen, 0.5);
  cone_radii_map["Uncorrected"]->SetMarkerStyle(21);
  cone_radii_map["Uncorrected"]->SetMarkerColor(kGreen);
  

  cone_radii_map["Corrected"]->Draw("E");
  cone_radii_map["Uncorrected"]->Draw("same E");  

  cone_radii_map["Corrected"]->GetXaxis()->SetTitle("pT");
  cone_radii_map["Uncorrected"]->GetXaxis()->CenterTitle();

  cone_radii_map["Corrected"]->GetYaxis()->SetRangeUser(10e-2, 10e6);
  cone_radii_map["Uncorrected"]->GetYaxis()->SetRangeUser(10e-2, 10e6);


  legend->AddEntry(cone_radii_map["Corrected"], "Corrected");
  legend->AddEntry(cone_radii_map["Uncorrected"], "Uncorrected");


  
  gPad->SetLogy();  

  legend->Draw();

  gPad->Print("ak5_distribution.pdf");

}

void invariant_mass() {

  double m_cut = 0.00;


  ifstream infile("../data/CMS_JetSample_invariant_mass.dat");

  TFile * rootFile_;
  
  string tag;
  int prescale;
  double uncorrected_invariant_mass, corrected_invariant_mass;

  unordered_map<std::string, TH1F * > invariant_mass;
  
  EColor colors[2] = {kRed, kGreen};
  
  gStyle->SetOptStat(false);


  invariant_mass["Corrected"] = new TH1F("", "", 50, 0.0, 50.0);
  invariant_mass["Uncorrected"] = new TH1F("", "", 50, 0.0, 50.0);


  string line;
  while(getline(infile, line)) {
    istringstream iss(line);
    
    iss >> tag >> uncorrected_invariant_mass >> corrected_invariant_mass >> prescale;
    
    if (tag != "#") {
      
      if (uncorrected_invariant_mass > m_cut)
        invariant_mass["Uncorrected"]->Fill(uncorrected_invariant_mass, prescale);

      if (corrected_invariant_mass > m_cut)
        invariant_mass["Corrected"]->Fill(corrected_invariant_mass, prescale);
    }
  }

  TCanvas *cst = new TCanvas("cst","Invariant Mass Spectrum", 1000, 600);

  TLegend * legend = new TLegend(0.7, 0.9, 0.9, 0.8);

  invariant_mass["Corrected"]->SetFillColorAlpha(kRed, 0.5);
  invariant_mass["Corrected"]->SetMarkerStyle(21);
  invariant_mass["Corrected"]->SetMarkerColor(kRed);
  

  invariant_mass["Uncorrected"]->SetFillColorAlpha(kGreen, 0.5);
  invariant_mass["Uncorrected"]->SetMarkerStyle(21);
  invariant_mass["Uncorrected"]->SetMarkerColor(kGreen);
  
  invariant_mass["Corrected"]->Draw("E");
  invariant_mass["Uncorrected"]->Draw("same E");  

  invariant_mass["Corrected"]->GetXaxis()->SetTitle("pT");
  invariant_mass["Uncorrected"]->GetXaxis()->CenterTitle();

  invariant_mass["Corrected"]->GetYaxis()->SetRangeUser(10e-2, 10e7);
  invariant_mass["Uncorrected"]->GetYaxis()->SetRangeUser(10e-2, 10e7);

  legend->AddEntry(invariant_mass["Corrected"], "Corrected");
  legend->AddEntry(invariant_mass["Uncorrected"], "Uncorrected");

  
  gPad->SetLogy();  

  legend->Draw();

  gPad->Print("invariant_mass_distribution.pdf");
}


void zg_plots() {

  double pt_cut = 153;


  ifstream infile("/media/aashish/opendata/eos/opendata/cms/Run2010B/Jet/analyzed.dat");

  TFile * rootFile_;
  TTree * multiplicityTree_;
  
  string tag, trigger_name;
  double uncorrected_pt, corrected_pt, zg_05, zg_1, zg_2;
  int prescale, event_number, run_number;

  unordered_map<double, TH1F * > z_cuts_map;
  vector<double> z_cuts = {0.05, 0.1, 0.2};

  EColor colors[3] = {kRed, kBlue, kGreen};
  const char * z_cuts_labels[3] = {"z_cut = 0.05", "z_cut = 0.1", "z_cut = 0.2"};
  
  gStyle->SetOptStat(false);

  for (unsigned int i = 0; i < z_cuts.size(); i++) {
      TH1F * cone_radius_temp = new TH1F("", "", 40, 0.0, 1.0);
      z_cuts_map[z_cuts[i]] = cone_radius_temp;
  }

  string line;
  while(getline(infile, line)) {
    istringstream iss(line);

    //# Entry        Event_Number          Run_Number          Hardest_pT     Corr_Hardest_pT            Prescale        Trigger_Name               zg_05                zg_1                zg_2
    
    iss >> tag >> event_number >> run_number >> uncorrected_pt >> corrected_pt >> prescale >> trigger_name >> zg_05 >> zg_1 >> zg_2;

    // cout << event_number << endl;
    // cout << zg_05 << zg_1 << zg_2 << prescale << endl;

    if (tag != "#") {
      if (corrected_pt >= pt_cut) {
          z_cuts_map[0.05]->Fill(zg_05, prescale);
          z_cuts_map[0.1]->Fill(zg_1, prescale);
          z_cuts_map[0.2]->Fill(zg_2, prescale);


      }
    }
  }

  TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

  for(unsigned int i = 0; i < z_cuts.size(); i++) {
      z_cuts_map[z_cuts[i]]->SetFillColorAlpha(colors[i], 0.5);
      z_cuts_map[z_cuts[i]]->SetMarkerStyle(21);
      z_cuts_map[z_cuts[i]]->SetMarkerColor(colors[i]);

      legend->AddEntry(z_cuts_map[z_cuts[i]], z_cuts_labels[i]);

      if (i == 0) {
        z_cuts_map[z_cuts[i]]->Draw("E");
        z_cuts_map[z_cuts[i]]->GetXaxis()->SetTitle("zg");
        z_cuts_map[z_cuts[i]]->GetXaxis()->CenterTitle();
        z_cuts_map[z_cuts[i]]->GetXaxis()->SetRangeUser(0., 0.6);
        z_cuts_map[z_cuts[i]]->GetYaxis()->SetRangeUser(0, 4000);
      }
      else {
        z_cuts_map[z_cuts[i]]->Draw("same E");            
      }
  }

  legend->Draw();

  gPad->Print("zg_distribution.pdf");
}


