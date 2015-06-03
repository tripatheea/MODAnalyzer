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

vector<string> split_string_to_components(string const &input);

void n_tilde_against_jet_multiplicity();
void hardest_pt_corresponding_triggers();
void fix_cone_radius_sweep_pt_cut();
void fix_pt_cut_sweep_cone_radius();

void plots() {
  fix_pt_cut_sweep_cone_radius();
}


void n_tilde_against_jet_multiplicity() {
  ifstream infile("../data/antikt_multiplicities.dat");

  TFile * rootFile_;
  TTree * multiplicityTree_;


  THStack *hs = new THStack("Fractional Jet Multiplicity", "Fractional Jet Multiplicity (pt_cut = 50.0 GeV, R = 0.5)");

  vector<TH1F * > N_tildes = vector<TH1F *>();
  EColor colors[7] = {kRed, kBlue, kGreen, kYellow, kMagenta, kOrange, kCyan};
  const char * antikT_labels[7] = {"antikT Jet Size = 0", "antikT Jet Size = 1", "antikT Jet Size = 2", "antikT Jet Size = 3", "antikT Jet Size = 4", "antikT Jet Size = 5", "antikT Jet Size = 6"};
  
  gStyle->SetOptStat(false);

  for (int i = 0; i < 7; i++) {
    TH1F * N_tilde_temp = new TH1F("", "", 50, -0.5, 6.0);
    N_tildes.push_back(N_tilde_temp);
  }
  
  double N_tilde;
  double antikt;

  double prescale_1, prescale_2;
  string name;

  string line;
  while(getline(infile, line)) {
    istringstream iss(line);
    vector<string> components = split_string_to_components(line);
    if (components[0] != "#") {
      double N_tilde = stod(components[2]);
      double jet_size = stod(components[3]);

      int prescale = stoi(components[6]) * stoi(components[7]);
      bool fired = (stoi(components[5]) == 1);
      
      double cone_radius = stod(components[8]);
      int pt_cut = stoi(components[9]);

      if ((fired) && (pt_cut == 50) && (cone_radius = 0.50)) {
        N_tildes[jet_size]->Fill(N_tilde, prescale);
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

  hs->Draw();
  hs->GetHistogram()->GetXaxis()->SetTitle("Fractional Jet Multiplicity (Jets Without Jets)");
  hs->GetHistogram()->GetXaxis()->CenterTitle();
  
  legend->Draw();
}



void hardest_pt_corresponding_triggers() {
  ifstream infile("../data/antikt_multiplicities.dat");

  TFile * rootFile_;
  TTree * multiplicityTree_;


  THStack *hs = new THStack("Hardest pt and corresponding trigger of jets", "Hardest pt and corresponding trigger of jets (pt_cut = 50.0 GeV, R = 0.5)");

  vector<TH1F * > hardest_pts = vector<TH1F *>();
  EColor colors[6] = {kRed, kBlue, kGreen, kYellow, kMagenta, kOrange};
  const char * trigger_labels[6] = {"HLT_Jet70U", "HLT_Jet50U", "HLT_Jet30U", "HLT_Jet15U", "HLT_L1Jet6U", "HLT_MinBiasPixel_SingleTrack"};
  
  gStyle->SetOptStat(false);

  for (int i = 0; i < 6; i++) {
    TH1F * pt_temp = new TH1F("a", "", 50, -0.5, 400);
    hardest_pts.push_back(pt_temp);
  }
  
  string line;
  while(getline(infile, line)) {
    istringstream iss(line);
    vector<string> components = split_string_to_components(line);
    if (components[0] != "#") {
      
      string trigger_name = components[4];

      int prescale = stoi(components[6]) * stoi(components[7]);
      bool fired = (stoi(components[5]) == 1);
      
      double cone_radius = stod(components[8]);
      int pt_cut = stoi(components[9]);
      double hardest_pt = stod(components[10]);

      if ((fired) && (pt_cut == 50) && (cone_radius = 0.50)) {
        int trigger_index = std::distance(trigger_labels, std::find(trigger_labels, trigger_labels + 6, trigger_name));
        hardest_pts[trigger_index]->Fill(hardest_pt, prescale);      
      }
    }
  }

  

  
  TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

  for(int i = 0; i < 6; i++) {
    hardest_pts[i]->SetFillColorAlpha(colors[i], 0.5);
    hardest_pts[i]->SetMarkerStyle(21);
    hardest_pts[i]->SetMarkerColor(colors[i]);
    hs->Add(hardest_pts[i]);
    legend->AddEntry(hardest_pts[i], trigger_labels[i]);
  }
  
  TCanvas *c2e = new TCanvas("c2e", "c2e", 600, 400);

  c2e->BuildLegend();

  gPad->SetLogy();

  hs->Draw();
  hs->GetHistogram()->GetXaxis()->SetTitle("pt_hardest");
  hs->GetHistogram()->GetXaxis()->CenterTitle();
  
  legend->Draw();
}




void fix_cone_radius_sweep_pt_cut() {

  // Fix R, sweep across pt_cut.

  double fixed_cone_radius = 0.5;

  ifstream infile("../data/antikt_multiplicities.dat");

  TFile * rootFile_;
  TTree * multiplicityTree_;


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
    vector<string> components = split_string_to_components(line);
    if (components[0] != "#") {
      
      double N_tilde = stod(components[2]);
      double jet_size = stod(components[3]);
      int prescale = stoi(components[6]) * stoi(components[7]);
      bool fired = (stoi(components[5]) == 1);
      
      double cone_radius = stod(components[8]);
      int pt_cut = stoi(components[9]);

      if ((fired) && (cone_radius == fixed_cone_radius)) {
        pt_cuts_map[pt_cut]->Fill(N_tilde, prescale);
      }
    }
  }


  TCanvas *cst = new TCanvas("cst","N_tilde and AntikT", 1000, 600);

  gPad->SetLogy();

  TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

  for(unsigned int i = 0; i < pt_cuts.size(); i++) {
    pt_cuts_map[pt_cuts[i]]->SetFillColorAlpha(colors[i], 0.5);
    pt_cuts_map[pt_cuts[i]]->SetMarkerStyle(21);
    pt_cuts_map[pt_cuts[i]]->SetMarkerColor(colors[i]);
    
    hs->Add(pt_cuts_map[pt_cuts[i]]);
  
    legend->AddEntry(pt_cuts_map[pt_cuts[i]], pt_cuts_labels[i]);
  }
  
  cst->BuildLegend();

  hs->Draw();

  hs->GetHistogram()->GetXaxis()->SetTitle("N_tilde");
  hs->GetHistogram()->GetXaxis()->CenterTitle();

  legend->Draw();

}





void fix_pt_cut_sweep_cone_radius() {

  // Fix pt_cut, sweep across R.
  int fixed_pt_cut = 80;

  ifstream infile("../data/antikt_multiplicities.dat");

  TFile * rootFile_;
  TTree * multiplicityTree_;


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
    vector<string> components = split_string_to_components(line);
    if (components[0] != "#") {
      
      double N_tilde = stod(components[2]);
      double jet_size = stod(components[3]);
      int prescale = stoi(components[6]) * stoi(components[7]);
      bool fired = (stoi(components[5]) == 1);
      
      double cone_radius = stod(components[8]);
      int pt_cut = stoi(components[9]);
      
      if ((fired) && (pt_cut == fixed_pt_cut)) {
        cone_radii_map[cone_radius]->Fill(N_tilde, prescale);
      }
    }
  }


  TCanvas *cst = new TCanvas("cst","N_tilde and AntikT", 1000, 600);

  gPad->SetLogy();

  TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

  for(unsigned int i = 0; i < cone_radii.size(); i++) {
    cone_radii_map[cone_radii[i]]->SetFillColorAlpha(colors[i], 0.5);
    cone_radii_map[cone_radii[i]]->SetMarkerStyle(21);
    cone_radii_map[cone_radii[i]]->SetMarkerColor(colors[i]);
    
    hs->Add(cone_radii_map[cone_radii[i]]);
  
    legend->AddEntry(cone_radii_map[cone_radii[i]], cone_radii_labels[i]);
  }

  cst->BuildLegend();

  hs->Draw();

  hs->GetHistogram()->GetXaxis()->SetTitle("N_tilde");
  hs->GetHistogram()->GetXaxis()->CenterTitle();

  legend->Draw();

}

vector<string> split_string_to_components(string const &input) { 
    istringstream buffer(input);
    vector<string> ret((istream_iterator<string>(buffer)), istream_iterator<string>());
    return ret;
}

