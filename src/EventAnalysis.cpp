R__LOAD_LIBRARY(TreeManager_C.so);

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "TreeManager.h"

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TMath.h>

// Function prototypes
TVector3 GetMPPC2DPos(int);
TVector3 GetMatchCubePos(int,int);
std::vector<TLorentzVector> Get3DMatchedCubes(std::vector<double>,double);
std::tuple<double,double,TGraph2D*> GetTrackLengthAngle(std::vector<TLorentzVector>);
bool CheckIndexMatch(int,int,int,int);
bool EventCosmicCut(std::vector<double>,double);

// Global parameters
double fCubeSize = 10; // In mm

void Event3DAnalysis(int file_option = 1, double ADC_cut = 500){

  std::string fin_name;

  if(file_option==1){
    fin_name = "../../inputs/GluedCubes_NoTeflon_NoGluedFiber.root";
  }
  else if(file_option==2){
    fin_name = "../../inputs/GluedCubes_WithTeflon_NoGluedFiber.root";
  }
  else if(file_option==3){
    fin_name = "../../inputs/GluedCubes_WithTeflon_GluedFiber.root";
  }

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();

  // Create some variables
  int n_pass = 0;
  bool pass_tag;
  std::vector<double> ADC_temp;
  std::vector<TLorentzVector> cube_array;
  double tot_ly;
  double trk_len;
  double trk_ang;
  TGraph2D *trk_graph;

  // Create some histograms
  // Total light yield vs. track direction
  TH2D *totly_ang = new TH2D("totly_ang","totly_ang",10,0,45,60,0,20000);
  totly_ang->GetXaxis()->SetTitle("Angle between track and Z axis / degree");
  totly_ang->GetYaxis()->SetTitle("Track light yield / ADC");
  totly_ang->GetXaxis()->SetLabelSize(0.04);
  totly_ang->GetXaxis()->SetTitleSize(0.04);
  totly_ang->GetYaxis()->SetLabelSize(0.04);
  totly_ang->GetYaxis()->SetTitleSize(0.04);
  totly_ang->GetYaxis()->SetTitleOffset(1.4);
  totly_ang->GetZaxis()->SetLabelSize(0.04);
  totly_ang->SetTitle("");

  // Distribution of local light yield per unit length 
  TH1D *locally = new TH1D("locally","locally",60,0,1000);
  locally->GetXaxis()->SetTitle("Average light yield ADC / mm");
  locally->GetYaxis()->SetTitle("Number of events / bin");
  locally->GetXaxis()->SetLabelSize(0.04);
  locally->GetXaxis()->SetTitleSize(0.04);
  locally->GetYaxis()->SetLabelSize(0.04);
  locally->GetYaxis()->SetTitleSize(0.04);
  locally->GetYaxis()->SetTitleOffset(1.4);
  locally->SetTitle("");

  // Loop over all events
  for(int n = 0; n < n_event; n++){

    data->GetMppc(n);

    // Check whether the event passed the cut
    ADC_temp.clear();
    for(int i = 0; i < 18; i++) ADC_temp.push_back(data->ADC(i));
    pass_tag = EventCosmicCut(ADC_temp,ADC_cut);

    if(pass_tag==false) continue;
    n_pass += 1;

    // Get the positions (+ light yields) of matched cubes
    // Each layer only consider the highest ADC in each plane, and use the center position of cube
    cube_array = Get3DMatchedCubes(ADC_temp,ADC_cut);

    // Get total light yield
    tot_ly = 0;
    for(int i = 0; i < cube_array.size(); i++) tot_ly += cube_array[i].T();

    // Get the total length of the track and angle with respect to positive Z axis
    std::tie(trk_len,trk_ang,trk_graph) = GetTrackLengthAngle(cube_array);

    // Fill the histograms
    totly_ang->Fill(trk_ang,tot_ly);
    locally->Fill(tot_ly/trk_len);

  }

  std::cout << "Total number of events: " << n_event << endl;
  std::cout << "Number of passed events: " << n_pass << endl;

  // Fit the local light yield distribution with Laudau function
  double range_low = locally->GetMean() - 2 * locally->GetRMS();
  double range_upp = locally->GetMean() + 4 * locally->GetRMS();
  locally->Fit("landau","","",range_low,range_upp);
  TF1 *fit_func = locally->GetFunction("landau");
  double mean = fit_func->GetParameter(1);
  double rms = fit_func->GetParameter(2);

  std::cout << "Local light yield mean: " << mean << ", RMS: " << rms << endl;

  TString name;
  TText *pl_name = new TText();
  pl_name->SetTextSize(0.04);
  TText *pl_mean = new TText();
  pl_mean->SetTextSize(0.03);
  TText *pl_rms = new TText();
  pl_rms->SetTextSize(0.03);
 
  // Draw the plots
  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("totly_ang","totly_ang",800,600);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  c1->cd();
  totly_ang->SetContour(99);
  totly_ang->Draw("colz");
  c1->Update();

  TCanvas *c2 = new TCanvas("locally","locally",700,600);
  c2->SetLeftMargin(0.15);
  c2->cd();
  locally->SetLineWidth(2);
  locally->SetLineColor(kBlue);
  locally->Draw("hist");

  name.Form("Landau fit:");
  pl_name->DrawTextNDC(0.55,0.68,name);
  name.Form("MPV = %f ADC / mm",mean);
  pl_mean->DrawTextNDC(0.55,0.61,name);
  //name.Form("Sigma = %f ADC / mm",rms);
  //pl_rms->DrawTextNDC(0.55,0.54,name);

  gPad->SetGridx();
  gPad->SetGridy();
  c2->Update();

  if(file_option==1){
    c1->SaveAs("../../../plots/scintillator_cube/gluedcubes_noteflon_nogluedfiber/totly_ang.png");
    c2->SaveAs("../../../plots/scintillator_cube/gluedcubes_noteflon_nogluedfiber/locally.png");
  }
  else if(file_option==2){
    c1->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_nogluedfiber/totly_ang.png");
    c2->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_nogluedfiber/locally.png");
  }
  else if(file_option==3){
    c1->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_gluedfiber/totly_ang.png");
    c2->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_gluedfiber/locally.png");
  }

  // Save the plots into output file
  TString fout_name;

  if(file_option==1){
    fout_name = "../../results/Event3DAnalysis_GluedCubes_NoTeflon_NoGluedFiber.root";
  }
  else if(file_option==2){
    fout_name = "../../results/Event3DAnalysis_GluedCubes_WithTeflon_NoGluedFiber.root";
  }
  else if(file_option==3){
    fout_name = "../../results/Event3DAnalysis_GluedCubes_WithTeflon_GluedFiber.root";
  }

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();

  c1->Write();
  c2->Write();

  totly_ang->Write();
  locally->Write();

  fout->Close();

}

// Draw event 3D display (after some cuts)
// Notation: file_option = input data file want to be used
// Notation: seed = random number generator
// Notation: ADC_cut = a cut applied on ADC value
void DrawEvent3D(int file_option = 1, int seed = 0, double ADC_cut = 500){

  std::string fin_name;

  if(file_option==1){
    fin_name = "../../inputs/GluedCubes_NoTeflon_NoGluedFiber.root";
  }
  else if(file_option==2){
    fin_name = "../../inputs/GluedCubes_WithTeflon_NoGluedFiber.root";
  }
  else if(file_option==3){
    fin_name = "../../inputs/GluedCubes_WithTeflon_GluedFiber.root";
  }

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();
  int rand_event;

  // Random number generator
  TRandom3 *rand = new TRandom3();
  rand->SetSeed(seed);

  // Create 3D plot to display event  
  TH3D *Event_show = new TH3D("Event_show","Event_show",3,0,3,3,0,3,3,0,3);
  Event_show->GetXaxis()->SetTitle("X axis cube");
  Event_show->GetYaxis()->SetTitle("Y axis cube");
  Event_show->GetZaxis()->SetTitle("Z axis cube");
  Event_show->SetTitle("");

  double size = 3 * fCubeSize;
  TH3D *bkg_temp = new TH3D("bkg_temp","bkg_temp",100,0,size,100,0,size,100,0,size);
  bkg_temp->GetXaxis()->SetTitle("X axis cube");
  bkg_temp->GetYaxis()->SetTitle("Y axis cube");
  bkg_temp->GetZaxis()->SetTitle("Z axis cube");
  bkg_temp->SetTitle("");

  // Create 2D plot to show MPPC on each plane
  TH2D *MPPC2D_xz = new TH2D("MPPC2D_xz","MPPC2D_xz",3,0,3,3,0,3);
  MPPC2D_xz->GetXaxis()->SetTitle("X axis cube");
  MPPC2D_xz->GetYaxis()->SetTitle("Z axis cube");
  MPPC2D_xz->GetXaxis()->SetTitleSize(0.04);
  MPPC2D_xz->GetXaxis()->SetLabelSize(0.04);
  MPPC2D_xz->GetYaxis()->SetTitleSize(0.04);
  MPPC2D_xz->GetYaxis()->SetLabelSize(0.04);
  MPPC2D_xz->GetYaxis()->SetTitleOffset(1.2);
  MPPC2D_xz->GetZaxis()->SetTitle("ADC");
  MPPC2D_xz->GetZaxis()->SetLabelSize(0.04);
  MPPC2D_xz->GetZaxis()->SetTitleOffset(1.2);
  MPPC2D_xz->SetTitle("MPPC XZ Projection");

  TH2D *MPPC2D_yz = (TH2D*)MPPC2D_xz->Clone("MPPC2D_yz");
  MPPC2D_yz->GetXaxis()->SetTitle("Y axis cube");
  MPPC2D_yz->SetTitle("MPPC YZ Projection");

  // Define some variables or vectors to save temporal data
  bool out_tag = false;
  bool pass_tag;
  double ADC_max;
  int index_xz, index_yz;
  double ADC_xz, ADC_yz;
  TVector3 cube_pos;
  std::vector<double> ADC_temp;
  std::vector<TLorentzVector> cube_array;
  TGraph2D *trk_graph;
  double trk_len; 
  double trk_ang;

  // Read the events to find one matches the cut
  while(out_tag==false){

    // Randomly choose an event
    rand_event = rand->Rndm() * n_event; 
    //std::cout << "Event number: " << rand_event << endl;
    data->GetMppc(rand_event);

    // Check whether this event passes the cosmic cut 
    ADC_temp.clear();
    for(int i = 0; i < 18; i++) ADC_temp.push_back(data->ADC(i));
    pass_tag = EventCosmicCut(ADC_temp,ADC_cut);

    if(pass_tag==true){

      // Get the cube array
      cube_array = Get3DMatchedCubes(ADC_temp,ADC_cut);

      // 3D plot show the cubes
      for(int i = 0; i < cube_array.size(); i++){
        Event_show->Fill(cube_array[i].X(),cube_array[i].Y(),cube_array[i].Z(),cube_array[i].T());
      }

      // 3D plot show the track (straight line)
      std::tie(trk_len,trk_ang,trk_graph) = GetTrackLengthAngle(cube_array);

      // Also draw 2D MPPC on each plane
      for(int i = 0; i < 18; i++){
        cube_pos = GetMPPC2DPos(i+1);
        MPPC2D_xz->Fill(cube_pos.X(),cube_pos.Z(),data->ADC(i));
        MPPC2D_yz->Fill(cube_pos.Y(),cube_pos.Z(),data->ADC(i)); 
      }
 
      out_tag = true; 
   
    }

  } 

  //trk_graph->GetXaxis()->SetRangeUser(0,3*fCubeSize);
  //trk_graph->GetYaxis()->SetRangeUser(0,3*fCubeSize);
  //trk_graph->GetZaxis()->SetRangeUser(0,3*fCubeSize);
  //trk_graph->GetXaxis()->SetTitle("X axis cube");
  //trk_graph->GetYaxis()->SetTitle("Y axis cube");
  //trk_graph->GetZaxis()->SetTitle("Z axis cube");
  //trk_graph->SetTitle("");
  trk_graph->SetLineWidth(2);
  trk_graph->SetLineColor(kBlue);

  // Plot the histograms
  gStyle->SetOptStat(0);
 
  TCanvas *c1 = new TCanvas("event_show","event_show",700,600);
  c1->cd();
  Event_show->Draw("BOX2");
  c1->Update();

  TCanvas *c2 = new TCanvas("MPPC2D_xz","MPPC2D_xz",700,600);
  c2->SetRightMargin(0.15);
  c2->cd();
  MPPC2D_xz->SetContour(99);
  MPPC2D_xz->Draw("colz TEXT45");
  c2->Update();

  TCanvas *c3 = new TCanvas("MPPC2D_yz","MPPC2D_yz",700,600);
  c3->SetRightMargin(0.15);
  c3->cd();
  MPPC2D_yz->SetContour(99);
  MPPC2D_yz->Draw("colz TEXT45");
  c3->Update();

  TCanvas *c4 = new TCanvas("combine","combine",1200,1200);
  c4->Divide(2,2);
  c4->cd(1);
  MPPC2D_xz->SetContour(99);
  MPPC2D_xz->Draw("colz TEXT45");
  c4->Update();
  c4->cd(2);
  MPPC2D_yz->SetContour(99);
  MPPC2D_yz->Draw("colz TEXT45");
  c4->Update();
  c4->cd(3);
  Event_show->Draw("BOX2");
  c4->cd(4);
  bkg_temp->Draw();
  trk_graph->Draw("LINE SAME");
  c4->Update();

  if(file_option==1){
    c4->SaveAs("../../../plots/scintillator_cube/gluedcubes_noteflon_nogluedfiber/combine.png");
  }
  else if(file_option==2){
    c4->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_nogluedfiber/combine.png");
  }
  else if(file_option==3){
    c4->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_gluedfiber/combine.png");
  }

  // Save the plot into output file
  TString fout_name;

  if(file_option==1){
    fout_name = "../../results/EventDisplay_GluedCubes_NoTeflon_NoGluedFiber.root";
  }
  else if(file_option==2){
    fout_name = "../../results/EventDisplay_GluedCubes_WithTeflon_NoGluedFiber.root";
  }
  else if(file_option==3){
    fout_name = "../../results/EventDisplay_GluedCubes_WithTeflon_GluedFiber.root";
  }

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();

  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();

  Event_show->Write();
  MPPC2D_xz->Write();
  MPPC2D_yz->Write();
  trk_graph->Write();

  fout->Close();

}

// Plot the ADC distribution for each channel
void DrawMPPCLightYield(int file_option = 1, double ADC_cut = 500){

  std::string fin_name;

  if(file_option==1){
    fin_name = "../../inputs/GluedCubes_NoTeflon_NoGluedFiber.root";
  }
  else if(file_option==2){
    fin_name = "../../inputs/GluedCubes_WithTeflon_NoGluedFiber.root";
  }
  else if(file_option==3){
    fin_name = "../../inputs/GluedCubes_WithTeflon_GluedFiber.root";
  }

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();

  // Create some variables
  std::vector<double> ADC_temp;
  bool passtag;

  // Light yield distributions for each channel 0 - 17
  TString name;
  TH1D *MPPC_ly[18];
  
  for(int i = 0; i < 18; i++){
    name.Form("MPPC_ly_chan%i",i+1);
    MPPC_ly[i] = new TH1D(name,name,80,500,4500);
    MPPC_ly[i]->GetXaxis()->SetTitle("ADC");
    MPPC_ly[i]->GetYaxis()->SetTitle("Number of events / bin");
    MPPC_ly[i]->SetTitle(name);
    MPPC_ly[i]->SetLineWidth(2);
    MPPC_ly[i]->SetLineColor(kBlue);
  }

  // Loop over all events
  for(int n = 0; n < n_event; n++){

    data->GetMppc(n);
 
    ADC_temp.clear();
    for(int i = 0; i < 18; i++) ADC_temp.push_back(data->ADC(i));
    passtag = EventCosmicCut(ADC_temp,ADC_cut);

    if(passtag==false) continue;  

    for(int i = 0; i < 18; i++) MPPC_ly[i]->Fill(data->ADC(i)); 

  }

  gStyle->SetOptStat(0);

  // MPPC channel face 1 + 3
  TCanvas *c1 = new TCanvas("MPPC2D_xz","MPPC2D_xz",1200,1200);
  c1->Divide(3,3);
  c1->cd(1);
  MPPC_ly[13]->Draw();
  //gPad->SetLogy();
  c1->cd(2);
  MPPC_ly[4]->Draw();
  //gPad->SetLogy();
  c1->cd(3);
  MPPC_ly[12]->Draw();
  //gPad->SetLogy();
  c1->cd(4);
  MPPC_ly[2]->Draw();
  //gPad->SetLogy(); 
  c1->cd(5);
  MPPC_ly[11]->Draw();
  //gPad->SetLogy();
  c1->cd(6);
  MPPC_ly[1]->Draw();
  //gPad->SetLogy();
  c1->cd(7);
  MPPC_ly[10]->Draw();
  //gPad->SetLogy();
  c1->cd(8);
  MPPC_ly[0]->Draw();
  //gPad->SetLogy();
  c1->cd(9);
  MPPC_ly[9]->Draw();
  //gPad->SetLogy();
  c1->Update();

  // MPPC channel face 2 + 4
  TCanvas *c2 = new TCanvas("MPPC2D_yz","MPPC2D_yz",1200,1200);
  c2->Divide(3,3);
  c2->cd(1);
  MPPC_ly[8]->Draw();
  //gPad->SetLogy();
  c2->cd(2);
  MPPC_ly[17]->Draw();
  //gPad->SetLogy();
  c2->cd(3);
  MPPC_ly[7]->Draw(); 
  //gPad->SetLogy();
  c2->cd(4);
  MPPC_ly[16]->Draw();
  //gPad->SetLogy();
  c2->cd(5);
  MPPC_ly[6]->Draw();
  //gPad->SetLogy();
  c2->cd(6);
  MPPC_ly[15]->Draw();
  //gPad->SetLogy();
  c2->cd(7);
  MPPC_ly[5]->Draw();
  //gPad->SetLogy();
  c2->cd(8);
  MPPC_ly[14]->Draw();
  //gPad->SetLogy();
  c2->cd(9);
  MPPC_ly[4]->Draw();
  //gPad->SetLogy();
  c2->Update();

  if(file_option==1){
    c1->SaveAs("../../../plots/scintillator_cube/gluedcubes_noteflon_nogluedfiber/MPPC2D_xz.png");
    c2->SaveAs("../../../plots/scintillator_cube/gluedcubes_noteflon_nogluedfiber/MPPC2D_yz.png");
  }
  else if(file_option==2){    
    c1->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_nogluedfiber/MPPC2D_xz.png");
    c2->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_nogluedfiber/MPPC2D_yz.png");
  }
  else if(file_option==3){
    c1->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_gluedfiber/MPPC2D_xz.png");
    c2->SaveAs("../../../plots/scintillator_cube/gluedcubes_withteflon_gluedfiber/MPPC2D_yz.png");
  }

  // Save the plots into output file
  TString fout_name;

  if(file_option==1){
    fout_name = "../../results/MPPCLightYield_GluedCubes_NoTeflon_NoGluedFiber.root";
  }
  else if(file_option==2){
    fout_name = "../../results/MPPCLightYield_GluedCubes_WithTeflon_NoGluedFiber.root";
  }
  else if(file_option==3){
    fout_name = "../../results/MPPCLightYield_GluedCubes_WithTeflon_GluedFiber.root";
  }

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();
  
  for(int i = 0; i < 18; i++) MPPC_ly[i]->Write();
 
  c1->Write();
  c2->Write();

  fout->Close();

} 

// ------------------------------------------------
// Below are auxiliary functions
// ------------------------------------------------

// Get the 3D matched cubes (one per layer)
// For each cube the 3D position and light yield will be returned
std::vector<TLorentzVector> Get3DMatchedCubes(std::vector<double> ADC_temp, double ADC_cut){

  double ADC_max;
  double ADC_xz, ADC_yz;
  int index_xz, index_yz;
  TVector3 cube_pos;
  TLorentzVector cube_temp;
  std::vector<TLorentzVector> cube_array;

  // In order to get 3D event, at each layer only consider the highest ADC on each plane
  
  // Top layer (XZ plane)
  ADC_max = 0;
  if(ADC_temp[13]>=ADC_max){
    ADC_max = ADC_temp[13];
    index_xz = 14;
    ADC_xz = ADC_temp[13];
  }
  if(ADC_temp[3]>=ADC_max){
    ADC_max = ADC_temp[3];
    index_xz = 4;
    ADC_xz = ADC_temp[3];
  }
  if(ADC_temp[12]>=ADC_max){
    ADC_max = ADC_temp[12];
    index_xz = 13;
    ADC_xz = ADC_temp[12];
  }
  // Top layer (YZ plane)
  ADC_max = 0;
  if(ADC_temp[8]>=ADC_max){
    ADC_max = ADC_temp[8];
    index_yz = 9;
    ADC_yz = ADC_temp[8];
  }
  if(ADC_temp[17]>=ADC_max){
    ADC_max = ADC_temp[17];
    index_yz = 18;
    ADC_yz = ADC_temp[17];
  }
  if(ADC_temp[7]>=ADC_max){
    ADC_max = ADC_temp[7];
    index_yz = 8;
    ADC_yz = ADC_temp[7];
  }
 
  cube_pos = GetMatchCubePos(index_xz,index_yz);
  cube_temp.SetXYZT(cube_pos.X(),cube_pos.Y(),cube_pos.Z(),ADC_xz+ADC_yz);
  cube_array.push_back(cube_temp);

  // Middle layer (XZ plane)
  ADC_max = 0;
  if(ADC_temp[2]>=ADC_max){
    ADC_max = ADC_temp[2];
    index_xz = 3;
    ADC_xz = ADC_temp[2];
  }
  if(ADC_temp[11]>=ADC_max){
    ADC_max = ADC_temp[11];
    index_xz = 12;
    ADC_xz = ADC_temp[11];
  }
  if(ADC_temp[1]>=ADC_max){
    ADC_max = ADC_temp[1];
    index_xz = 2;
    ADC_xz = ADC_temp[1];
  }
  // Middle layer (YZ plane)
  ADC_max = 0;
  if(ADC_temp[16]>=ADC_max){
    ADC_max = ADC_temp[16];
    index_yz = 17;
    ADC_yz = ADC_temp[16];
  }
  if(ADC_temp[6]>=ADC_max){
    ADC_max = ADC_temp[6];
    index_yz = 7;
    ADC_yz = ADC_temp[6];
  }
  if(ADC_temp[15]>=ADC_max){
    ADC_max = ADC_temp[15];
    index_yz = 16;
    ADC_yz = ADC_temp[15];
  }

  cube_pos = GetMatchCubePos(index_xz,index_yz);
  cube_temp.SetXYZT(cube_pos.X(),cube_pos.Y(),cube_pos.Z(),ADC_xz+ADC_yz);
  cube_array.push_back(cube_temp);

  // Bottom layer (XZ plane)
  ADC_max = 0;
  if(ADC_temp[10]>=ADC_max){
    ADC_max = ADC_temp[10];
    index_xz = 11;
    ADC_xz = ADC_temp[10];
  }
  if(ADC_temp[0]>=ADC_max){
    ADC_max = ADC_temp[0];
    index_xz = 1;
    ADC_xz = ADC_temp[0];
  }
  if(ADC_temp[9]>=ADC_max){
    ADC_max = ADC_temp[9];
    index_xz = 10;
    ADC_xz = ADC_temp[9];
  }

  // Bottom layer (YZ plane)
  ADC_max = 0;
  if(ADC_temp[5]>=ADC_max){
    ADC_max = ADC_temp[5];
    index_yz = 6;
    ADC_yz = ADC_temp[5];
  }
  if(ADC_temp[14]>=ADC_max){
    ADC_max = ADC_temp[14];
    index_yz = 15;
    ADC_yz = ADC_temp[14];
  }
  if(ADC_temp[4]>=ADC_max){
    ADC_max = ADC_temp[4];
    index_yz = 5;
    ADC_yz = ADC_temp[4];
  }

  cube_pos = GetMatchCubePos(index_xz,index_yz);
  cube_temp.SetXYZT(cube_pos.X(),cube_pos.Y(),cube_pos.Z(),ADC_xz+ADC_yz);
  cube_array.push_back(cube_temp);

  return cube_array;

}

// Estimate the track length in 3*3*3 cube matrix and the angle with respect to positive Z axis
std::tuple<double,double,TGraph2D*> GetTrackLengthAngle(std::vector<TLorentzVector> cube_array){

  // Create some parameters
  TGraph *fit_xz = new TGraph();
  TGraph *fit_yz = new TGraph();
  TF1 *fit_func = new TF1("fit_func","[0]*x+[1]");
  TF1 *fit_back;
  TGraph2D *graph = new TGraph2D();
  double x_temp, y_temp, z_temp;
  double ax, bx;
  double ay, by;
  double dx, dy;
  double d;
  double ang;
  double bias = fCubeSize / 2;
  double size = 3 * fCubeSize;

  // Fill the graph
  for(int i = 0; i < cube_array.size(); i++){
 
    // XZ plane
    x_temp = cube_array[i].Z() * fCubeSize + bias;
    y_temp = cube_array[i].X() * fCubeSize + bias;
    fit_xz->SetPoint(i,x_temp,y_temp);

    // YZ plane
    x_temp = cube_array[i].Z() * fCubeSize + bias;
    y_temp = cube_array[i].Y() * fCubeSize + bias;
    fit_yz->SetPoint(i,x_temp,y_temp);

    //std::cout << cube_array[i].X() << " " << cube_array[i].Y() << " " << cube_array[i].Z() << endl;

  }

  // Fit the graph with straight line
  fit_xz->Fit("fit_func","","",0,size);
  fit_back = fit_xz->GetFunction("fit_func");
  ax = fit_back->GetParameter(0);
  bx = fit_back->GetParameter(1);
  dx = size * ax;

  fit_yz->Fit("fit_func","","",0,size);
  fit_back = fit_yz->GetFunction("fit_func");
  ay = fit_back->GetParameter(0);
  by = fit_back->GetParameter(1);
  dy = size * ay;
 
  // Compute the track total length
  d = sqrt(pow(dx,2) + pow(dy,2) + pow(size,2));

  // Compute the track angle with respect to positive Z axis
  ang = TMath::ACos(size/d) * 180 / TMath::Pi();

  // Draw the track (straight line)
  const int npoint = 100;
  double step = size / npoint;

  for(int i = 0; i < npoint; i++){
    z_temp = step * i;
    x_temp = ax * z_temp + bx;
    y_temp = ay * z_temp + by;
    graph->SetPoint(i,x_temp,y_temp,z_temp);
  }

  return std::make_tuple(d,ang,graph);

}

// Check whether the event passes the cosmic cut
bool EventCosmicCut(std::vector<double> ADC_temp, double ADC_cut){

  // Define some variables  
  int n_topxzhit = 0;
  //int n_topxzhit_cut = 0;
  int n_topyzhit = 0;
  //int n_topyzhit_cut = 0;
  int n_botxzhit = 0;
  //int n_botxzhit_cut = 0;
  int n_botyzhit = 0;
  //int n_botyzhit_cut = 0;
  int n_tothit = 0;
  
  bool isHit; 

  // Loop over all MPPC channels (currently 18)
  for(int i = 1; i <= 18; i++){

    isHit = false;

    // Check if ADC higher than the cut
    if(ADC_temp[i-1]>=ADC_cut){
      isHit = true;
      n_tothit += 1;
    } 

    if(isHit==false) continue;

    // Top layer (XZ plane)
    if(i==14 || i==4 || i==13) n_topxzhit += 1;
    // Top layer (YZ plane)
    else if(i==9 || i==18 || i==8) n_topyzhit += 1;
    // Bottom layer (XZ plane)
    else if(i==11 || i==1 || i==10) n_botxzhit += 1;
    // Bottom layer (YZ plane)
    else if(i==6 || i==15 || i==5) n_botyzhit += 1;
      
  }

  // Check if the requirements are satisfied 
  if(n_tothit>=6 && n_topxzhit==1 && n_topyzhit==1 && n_botxzhit==1 && n_botyzhit==1){
    return true;
  }
  else{
    return false;
  }

}

// The test cubes are a 3*3*3 matrix
TVector3 GetMPPC2DPos(int MPPCChan){

  TVector3 pos;

  // Face 1
  if(MPPCChan==1) pos.SetXYZ(1,-1,0);
  else if(MPPCChan==2) pos.SetXYZ(2,-1,1);
  else if(MPPCChan==3) pos.SetXYZ(0,-1,1);
  else if(MPPCChan==4) pos.SetXYZ(1,-1,2);
  // Face 3
  else if(MPPCChan==10) pos.SetXYZ(0,-1,0);
  else if(MPPCChan==11) pos.SetXYZ(2,-1,0);
  else if(MPPCChan==12) pos.SetXYZ(1,-1,1);
  else if(MPPCChan==13) pos.SetXYZ(0,-1,2);
  else if(MPPCChan==14) pos.SetXYZ(2,-1,2);
  // Face 2
  else if(MPPCChan==5) pos.SetXYZ(-1,2,0);
  else if(MPPCChan==6) pos.SetXYZ(-1,0,0);
  else if(MPPCChan==7) pos.SetXYZ(-1,1,1);
  else if(MPPCChan==8) pos.SetXYZ(-1,2,2);
  else if(MPPCChan==9) pos.SetXYZ(-1,0,2);
  // Face 4
  else if(MPPCChan==15) pos.SetXYZ(-1,1,0);
  else if(MPPCChan==16) pos.SetXYZ(-1,0,1);
  else if(MPPCChan==17) pos.SetXYZ(-1,2,1);
  else if(MPPCChan==18) pos.SetXYZ(-1,1,2); 
  // No position available
  else pos.SetXYZ(-1,-1,-1);

  return pos;

}

// Return the cube position where two fibers x and y intersect
TVector3 GetMatchCubePos(int x, int y){

  TVector3 pos;

  // Bottom layer
  if(CheckIndexMatch(x,y,6,10)==true) pos.SetXYZ(0,0,0);
  else if(CheckIndexMatch(x,y,6,1)==true) pos.SetXYZ(1,0,0);
  else if(CheckIndexMatch(x,y,6,11)==true) pos.SetXYZ(2,0,0);
  else if(CheckIndexMatch(x,y,15,10)==true) pos.SetXYZ(0,1,0);
  else if(CheckIndexMatch(x,y,15,1)==true) pos.SetXYZ(1,1,0);
  else if(CheckIndexMatch(x,y,15,11)==true) pos.SetXYZ(2,1,0);
  else if(CheckIndexMatch(x,y,5,10)==true) pos.SetXYZ(0,2,0);
  else if(CheckIndexMatch(x,y,5,1)==true) pos.SetXYZ(1,2,0);
  else if(CheckIndexMatch(x,y,5,11)==true) pos.SetXYZ(2,2,0);
  // Middle layer
  else if(CheckIndexMatch(x,y,16,3)==true) pos.SetXYZ(0,0,1);
  else if(CheckIndexMatch(x,y,16,12)==true) pos.SetXYZ(1,0,1);
  else if(CheckIndexMatch(x,y,16,2)==true) pos.SetXYZ(2,0,1);
  else if(CheckIndexMatch(x,y,7,3)==true) pos.SetXYZ(0,1,1);
  else if(CheckIndexMatch(x,y,7,12)==true) pos.SetXYZ(1,1,1);
  else if(CheckIndexMatch(x,y,7,2)==true) pos.SetXYZ(2,1,1);
  else if(CheckIndexMatch(x,y,17,3)==true) pos.SetXYZ(0,2,1);
  else if(CheckIndexMatch(x,y,17,12)==true) pos.SetXYZ(1,2,1);
  else if(CheckIndexMatch(x,y,17,2)==true) pos.SetXYZ(2,2,1);
  // Top layer
  else if(CheckIndexMatch(x,y,9,13)==true) pos.SetXYZ(0,0,2);
  else if(CheckIndexMatch(x,y,9,4)==true) pos.SetXYZ(1,0,2);
  else if(CheckIndexMatch(x,y,9,14)==true) pos.SetXYZ(2,0,2);
  else if(CheckIndexMatch(x,y,18,13)==true) pos.SetXYZ(0,1,2);
  else if(CheckIndexMatch(x,y,18,4)==true) pos.SetXYZ(1,1,2);
  else if(CheckIndexMatch(x,y,18,14)==true) pos.SetXYZ(2,1,2);
  else if(CheckIndexMatch(x,y,8,13)==true) pos.SetXYZ(0,2,2);
  else if(CheckIndexMatch(x,y,8,4)==true) pos.SetXYZ(1,2,2);
  else if(CheckIndexMatch(x,y,8,14)==true) pos.SetXYZ(2,2,2);
  // Not matched
  else pos.SetXYZ(-1,-1,-1);
  
  return pos;

}

bool CheckIndexMatch(int x, int y, int a, int b){

  if(x==a && y==b) return true;
  else if(x==b && y==a) return true;
  else return false;

}
