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


// void fractional_multiplicities() {
//   std::ifstream infile("antikt_multiplicities.csv");

//   TFile * rootFile_;
//   TTree * multiplicityTree_;

//   TCanvas *c2e = new TCanvas("ntilde_antikt", "Ntilde and Antikt Multiplicity", 1000, 600);
  
//   gStyle->SetOptStat(false);
//   TH1F *h1 = new TH1F("Antikt Multiplicity", "", 50, -0.5, 6.0);
//   TH1F *h2 = new TH1F("Ntilde", "", 50, -0.5, 6.0);

//   double N_tilde;
//   double antikt;
//   double prescale_1, prescale_2;

//   while(infile >> N_tilde >> antikt >> prescale_1 >> prescale_2) {
//     h1->Fill(antikt, prescale_1 * prescale_2);
//     h2->Fill(N_tilde, prescale_1 * prescale_2);
//   }


//   h1->Draw();
//   c2e->Update();

//   // Scale h2 to the pad coordinates.
//   Float_t rightmax = 1.1 * h2->GetMaximum();
//   Float_t scale = gPad->GetUymax() / rightmax;

//   h2->SetLineColor(kRed);
//   h2->Scale(scale);
//   h2->Draw("same");

//   gPad->SetLogy();

//   // Draw an axis on the right side.
//   TGaxis *axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), gPad->GetUymax(), 0, rightmax, 510, "+L");
//   axis->SetLineColor(kRed);
//   axis->SetLabelColor(kRed);
//   axis->Draw();

//   c2e->Update();

//   TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
//   leg->AddEntry(h1, "Antikt");
//   leg->AddEntry(h2, "Ntilde");
//   leg->Draw();

// }



// void fractional_multiplicities() {
//   std::ifstream infile("antikt_multiplicities.csv");

//   TFile * rootFile_;
//   TTree * multiplicityTree_;

//   // Create, fill and project a 2D histogram.
//   TCanvas *c2e = new TCanvas("c2e", "c2e", 600, 400);
//   TH2F *h2 = new TH2F("Antikt Multiplicity and Ntilde", "", 50, -0.5, 6.0, 50, -0.5, 6.0);

//   rootFile_ = new TFile("antikt_multiplicities.root", "RECREATE");
//   multiplicityTree_ = new TTree("Multiplicity", "Multiplicity");
  
//   double N_tilde;
//   int antikt;

//   double prescale_1, prescale_2;

//   multiplicityTree_->Branch("multiplicity", &N_tilde, "N_tilde/D");
//   multiplicityTree_->Branch("antikt", &antikt, "antikt/I");

//   while(infile >> N_tilde >> antikt >> prescale_1 >> prescale_2) {
//     multiplicityTree_->Fill();
//     h2->Fill(N_tilde, antikt, prescale_1 * prescale_2);
//   }

//   rootFile_->cd();
  
//   multiplicityTree_->Write();

//   rootFile_->Close();



//   // Draw.
//   TH1D * projh2X = h2->ProjectionX();
//   TH1D * projh2Y = h2->ProjectionY();
  
//   h2->Draw();

// }



// void fractional_multiplicities() {
//   std::ifstream infile("antikt_multiplicities.csv");

//   TFile * rootFile_;
//   TTree * multiplicityTree_;

//   TH1F *h1st = new TH1F("h1st","antikt",50, -0.25, 4);
//   TH1F *h2st = new TH1F("h2st","N_tilde",50, -0.25, 4); 

//    double N_tilde;
//   double antikt;
//   double prescale_1, prescale_2;
  
//   while(infile >> N_tilde >> antikt >> prescale_1 >> prescale_2) {
//     h1st->Fill(antikt, prescale_1 * prescale_2);
//     h2st->Fill(N_tilde, prescale_1 * prescale_2);
//   }

//   THStack *hs = new THStack("ntilde_antikt","Ntilde and AntikT Multiplicity");



//   h1st->SetFillColor(kRed);
//   h1st->SetMarkerStyle(21);
//   h1st->SetMarkerColor(kRed);
//   hs->Add(h1st);

//   h2st->SetFillColor(kBlue);
//   h2st->SetMarkerStyle(21);
//   h2st->SetMarkerColor(kBlue);
//   hs->Add(h2st);


//   gStyle->SetOptStat(false);
//   TCanvas *cst = new TCanvas("cst","N_tilde and AntikT", 1000, 600);

//     gPad->SetLogy();

//   hs->Draw();

//   TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
//   leg->AddEntry(h1st, "Antikt");
//   leg->AddEntry(h2st, "N_tilde");
//   leg->Draw();  
// }



// void fractional_multiplicities() {
//   std::ifstream infile("antikt_multiplicities.csv");

//   TFile * rootFile_;
//   TTree * multiplicityTree_;


//   THStack *hs = new THStack("Fractional Jet Multiplicity", "Fractional Jet Multiplicity (pt_cut = 50.0 GeV, R = 0.5)");

//   vector<TH1F * > N_tildes = vector<TH1F *>();
//   EColor colors[7] = {kRed, kBlue, kGreen, kYellow, kMagenta, kOrange, kCyan};
//   const char * antikT_labels[7] = {"antikT Jet Size = 0", "antikT Jet Size = 1", "antikT Jet Size = 2", "antikT Jet Size = 3", "antikT Jet Size = 4", "antikT Jet Size = 5", "antikT Jet Size = 6"};
  
//   gStyle->SetOptStat(false);

//   for (int i = 0; i < 7; i++) {
//     TH1F * N_tilde_temp = new TH1F("", "", 50, -0.5, 6.0);
//     N_tildes.push_back(N_tilde_temp);
//   }
  
//   double N_tilde;
//   double antikt;

//   double prescale_1, prescale_2;
//   string name;

//   while(infile >> N_tilde >> antikt >> prescale_1 >> prescale_2 >> name) {
//     cout << N_tilde << endl;
//     N_tildes[antikt]->Fill(N_tilde);
//   }

//   TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

//   for(int i = 0; i < 6; i++) {
//     N_tildes[i]->SetFillColorAlpha(colors[i], 0.5);
//     N_tildes[i]->SetMarkerStyle(21);
//     N_tildes[i]->SetMarkerColor(colors[i]);
//     hs->Add(N_tildes[i]);
//     legend->AddEntry(N_tildes[i], antikT_labels[i]);
//   }
  
//   TCanvas *c2e = new TCanvas("c2e", "c2e", 600, 400);

//   c2e->BuildLegend();

//   gPad->SetLogy();

//   hs->Draw();
//   hs->GetHistogram()->GetXaxis()->SetTitle("Fractional Jet Multiplicity (Jets Without Jets)");
//   hs->GetHistogram()->GetXaxis()->CenterTitle();
  
//   legend->Draw();
// }



void fractional_multiplicities() {
  ifstream infile("antikt_multiplicities.dat");

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
      double hardest_pt = stod(components[10]);
      bool fired = (stoi(components[5]) == 1);

      


      if (fired) {
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




// void fractional_multiplicities() {

//   // Fix R, sweep across pt_cut.

//   ifstream infile("antikt_multiplicities.dat");

//   TFile * rootFile_;
//   TTree * multiplicityTree_;


//   THStack *hs = new THStack("Fractional Jet Multiplicity", "Fractional Jet Multiplicity (R = 0.5)");



//   unordered_map<double, TH1F * > pt_cuts_map;
//   vector<int> pt_cuts = {50, 80, 110};

//   EColor colors[3] = {kRed, kBlue, kGreen};
//   const char * pt_cuts_labels[3] = {"pT_cut = 50 GeV", "pT_cut = 80 GeV", "pT_cut = 110 GeV"};
  
//   gStyle->SetOptStat(false);

//   for (unsigned int i = 0; i < pt_cuts.size(); i++) {
//     TH1F * pt_cut_temp = new TH1F("", "", 50, -0.5, 6.0);
//     pt_cuts_map[pt_cuts[i]] = pt_cut_temp;
//   }


//   string line;
//   while(getline(infile, line)) {
//     istringstream iss(line);
//     vector<string> components = split_string_to_components(line);
//     if (components[0] != "#") {
      
//       double N_tilde = stod(components[2]);
//       double jet_size = stod(components[3]);
//       int prescale = stoi(components[6]) * stoi(components[7]);
//       bool fired = (stoi(components[5]) == 1);
      
//       double cone_radius = stod(components[8]);
//       int pt_cut = stoi(components[9]);
      
//       if (fired) {
//         pt_cuts_map[pt_cut]->Fill(N_tilde, prescale);
//       }
//     }
//   }


//   TCanvas *cst = new TCanvas("cst","N_tilde and AntikT", 1000, 600);

//   // gPad->SetLogy();

//   TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

//   for(unsigned int i = 0; i < pt_cuts.size(); i++) {
//     pt_cuts_map[pt_cuts[i]]->SetFillColorAlpha(colors[i], 0.5);
//     pt_cuts_map[pt_cuts[i]]->SetMarkerStyle(21);
//     pt_cuts_map[pt_cuts[i]]->SetMarkerColor(colors[i]);
    
//     hs->Add(pt_cuts_map[pt_cuts[i]]);
  
//     legend->AddEntry(pt_cuts_map[pt_cuts[i]], pt_cuts_labels[i]);
//   }
  


//   cst->BuildLegend();

//   hs->Draw();

//   hs->GetHistogram()->GetXaxis()->SetTitle("N_tilde");
//   hs->GetHistogram()->GetXaxis()->CenterTitle();

//   legend->Draw();

// }





// void fractional_multiplicities() {

//   // Fix pt_cut, sweep across R.

//   ifstream infile("antikt_multiplicities.dat");

//   TFile * rootFile_;
//   TTree * multiplicityTree_;


//   THStack *hs = new THStack("Fractional Jet Multiplicity", "Fractional Jet Multiplicity (pt_cut = 50 GeV)");



//   unordered_map<double, TH1F * > cone_radii_map;
//   vector<double> cone_radii = {0.3, 0.5, 0.7};

//   EColor colors[3] = {kRed, kBlue, kGreen};
//   const char * cone_radii_labels[3] = {"R = 0.3", "R = 0.5", "R = 0.7"};
  
//   gStyle->SetOptStat(false);

//   for (unsigned int i = 0; i < cone_radii.size(); i++) {
//     TH1F * cone_radius_temp = new TH1F("", "", 50, -0.5, 6.0);
//     cone_radii_map[cone_radii[i]] = cone_radius_temp;
//   }


//   string line;
//   while(getline(infile, line)) {
//     istringstream iss(line);
//     vector<string> components = split_string_to_components(line);
//     if (components[0] != "#") {
      
//       double N_tilde = stod(components[2]);
//       double jet_size = stod(components[3]);
//       int prescale = stoi(components[6]) * stoi(components[7]);
//       bool fired = (stoi(components[5]) == 1);
      
//       double cone_radius = stod(components[8]);
//       int pt_cut = stoi(components[9]);
      
//       if (fired) {
//         cone_radii_map[cone_radius]->Fill(N_tilde, prescale);
//       }
//     }
//   }


//   TCanvas *cst = new TCanvas("cst","N_tilde and AntikT", 1000, 600);

//   // gPad->SetLogy();

//   TLegend * legend = new TLegend(0.6, 0.7, 0.85, 0.9);

//   for(unsigned int i = 0; i < cone_radii.size(); i++) {
//     cone_radii_map[cone_radii[i]]->SetFillColorAlpha(colors[i], 0.5);
//     cone_radii_map[cone_radii[i]]->SetMarkerStyle(21);
//     cone_radii_map[cone_radii[i]]->SetMarkerColor(colors[i]);
    
//     hs->Add(cone_radii_map[cone_radii[i]]);
  
//     legend->AddEntry(cone_radii_map[cone_radii[i]], cone_radii_labels[i]);
//   }
  


//   cst->BuildLegend();

//   hs->Draw();

//   hs->GetHistogram()->GetXaxis()->SetTitle("N_tilde");
//   hs->GetHistogram()->GetXaxis()->CenterTitle();

//   legend->Draw();

// }


vector<string> split_string_to_components(string const &input) { 
    istringstream buffer(input);
    vector<string> ret((istream_iterator<string>(buffer)), istream_iterator<string>());
    return ret;
}