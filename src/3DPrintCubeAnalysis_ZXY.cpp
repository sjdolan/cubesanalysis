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
#include <TGraph2D.h>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TMath.h>

// Function prototypes
TLorentzVector GetNearbyChannel(int);
std::tuple<bool,int> CheckEventCosmicVertical(Mppc*);
TVector3 GetTopChanPos(int);
TVector3 GetTrigChanPos(int,int);
int ReturnIndex(int);

// Global parameters
std::string fin_name = "../../inputs/3dprinted_cubes/3DPrintCubes_NewBoard.root";

const int fChanNum = 9; // 3*3 matrix
int fChanOrder[fChanNum] = {15,17,19,21,23,25,27,29,31}; // 3*3 matrix channel number

double f13TopChanOrder[3] = {4,2,6}; // Channel number for trigger cubes
double f13DownChanOrder[3] = {3,1,5};
double f24TopChanOrder[3] = {12,14,8};
double f24DownChanOrder[3] = {11,13,7};

const double fADCCut = 500;

// ---------------------------
// Main functions
// ---------------------------

void CrosstalkAnalysis(){

  bool swap = false;
  bool newswap = false;

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();
  
  std::cout << "Total number of events: " << n_event << std::endl;
  
  TString name;
  TH1D *noise_esti[fChanNum];
  TH1D *xtalk_esti[fChanNum];
  TH1D *xtalk_adc[fChanNum];
  TH1D *track_adc[fChanNum];
  
  for(int i = 0; i < fChanNum; i++){
    name.Form("channel%i_noise",fChanOrder[i]);
    noise_esti[i] = new TH1D(name,name,100,1,400);
    noise_esti[i]->GetXaxis()->SetTitle("Estimated noise level / ADC");
    noise_esti[i]->GetYaxis()->SetTitle("Number of events / bin");
    noise_esti[i]->GetXaxis()->SetLabelSize(0.04);
    noise_esti[i]->GetXaxis()->SetTitleSize(0.04);
    noise_esti[i]->GetYaxis()->SetLabelSize(0.04);
    noise_esti[i]->GetYaxis()->SetTitleSize(0.04);
    noise_esti[i]->GetYaxis()->SetTitleOffset(1.4);
    noise_esti[i]->SetTitle(name);
    noise_esti[i]->SetLineWidth(2);
    noise_esti[i]->SetLineColor(kBlue);
    
    name.Form("channel%i_xtalk",fChanOrder[i]);
    xtalk_esti[i] = new TH1D(name,name,30,0,15);
    xtalk_esti[i]->GetXaxis()->SetTitle("Crosstalk fraction / %");
    xtalk_esti[i]->GetYaxis()->SetTitle("Number of events / bin");
    xtalk_esti[i]->GetXaxis()->SetLabelSize(0.04);
    xtalk_esti[i]->GetXaxis()->SetTitleSize(0.04);
    xtalk_esti[i]->GetYaxis()->SetLabelSize(0.04);
    xtalk_esti[i]->GetYaxis()->SetTitleSize(0.04);
    xtalk_esti[i]->GetYaxis()->SetTitleOffset(1.4);
    xtalk_esti[i]->SetTitle(name);
    xtalk_esti[i]->SetLineWidth(2);
    xtalk_esti[i]->SetLineColor(kBlue);
    
    name.Form("channel%i_xtalkADC",fChanOrder[i]);
    xtalk_adc[i] = new TH1D(name,name,100,0,1000);
    xtalk_adc[i]->GetXaxis()->SetTitle("Crosstalk channel / ADC");
    xtalk_adc[i]->GetYaxis()->SetTitle("Number of events / bin");
    xtalk_adc[i]->GetXaxis()->SetLabelSize(0.04);
    xtalk_adc[i]->GetXaxis()->SetTitleSize(0.04);
    xtalk_adc[i]->GetYaxis()->SetLabelSize(0.04);
    xtalk_adc[i]->GetYaxis()->SetTitleSize(0.04);
    xtalk_adc[i]->GetYaxis()->SetTitleOffset(1.4);
    xtalk_adc[i]->SetTitle(name);
    xtalk_adc[i]->SetLineWidth(2);
    xtalk_adc[i]->SetLineColor(kBlue);

    name.Form("channel%i_trackADC",fChanOrder[i]);
    track_adc[i] = new TH1D(name,name,100,0,4000);
    track_adc[i]->GetXaxis()->SetTitle("Track channel / ADC");
    track_adc[i]->GetYaxis()->SetTitle("Number of events / bin");
    track_adc[i]->GetXaxis()->SetLabelSize(0.04);
    track_adc[i]->GetXaxis()->SetTitleSize(0.04);
    track_adc[i]->GetYaxis()->SetLabelSize(0.04);
    track_adc[i]->GetYaxis()->SetTitleSize(0.04);
    track_adc[i]->GetYaxis()->SetTitleOffset(1.4);
    track_adc[i]->SetTitle(name);
    track_adc[i]->SetLineWidth(2);
    track_adc[i]->SetLineColor(kBlue);
  }
 
  TH1D *xtalk_all = new TH1D("xtalk_rate","xtalk_rate",30,0,15);
  xtalk_all->GetXaxis()->SetTitle("Crosstalk fraction / %");
  xtalk_all->GetYaxis()->SetTitle("Number of events / bin");
  xtalk_all->GetXaxis()->SetLabelSize(0.05);
  xtalk_all->GetXaxis()->SetTitleSize(0.05);
  xtalk_all->GetYaxis()->SetLabelSize(0.05);
  xtalk_all->GetYaxis()->SetTitleSize(0.05);
  xtalk_all->GetYaxis()->SetTitleOffset(1.4);
  xtalk_all->SetTitle("");
  
  // Estimate noise level
  std::cout << "Start to estimate noise level" << std::endl;
  for(int n = 0; n < n_event; n++){

    data->GetMppc(n,swap,newswap);
    
    for(int i = 0; i < fChanNum; i++){
      double cen_adc = data->ADC(fChanOrder[i]-1);
      if(cen_adc>fADCCut) continue;
      
      TLorentzVector chan_near = GetNearbyChannel(fChanOrder[i]);
      
      bool nearby_cut = true;
      if(chan_near.X()!=-1){
        if(data->ADC(chan_near.X()-1)>fADCCut) nearby_cut = false;
      }
      if(chan_near.Y()!=-1){
        if(data->ADC(chan_near.Y()-1)>fADCCut) nearby_cut = false;
      }
      if(chan_near.Z()!=-1){
        if(data->ADC(chan_near.Z()-1)>fADCCut) nearby_cut = false;
      }
      if(chan_near.T()!=-1){
        if(data->ADC(chan_near.T()-1)>fADCCut) nearby_cut = false;
      }
      
      if(nearby_cut==false) continue;
      
      noise_esti[i]->Fill(cen_adc);
    }

  }
  
  double noise_mean[fChanNum];
  double noise_rms[fChanNum];
  for(int i = 0; i < fChanNum; i++){
    noise_mean[i] = noise_esti[i]->GetMean();
    noise_rms[i] = noise_esti[i]->GetRMS();
  }
  
  // Then estimate crosstalk
  for(int n = 0; n < n_event; n++){
  
    data->GetMppc(n,swap,newswap);
    
    // Only select vertical tracks and should pass cosmic cut
    bool pass_tag;
    int trk_index;
    std::tie(pass_tag,trk_index) = CheckEventCosmicVertical(data);
  
    if(pass_tag!=true) continue;

    double cen_adc = data->ADC(fChanOrder[trk_index]-1);
    track_adc[trk_index]->Fill(cen_adc);
    
    TLorentzVector near_chan = GetNearbyChannel(fChanOrder[trk_index]);
    
    double xtalk_rate;
    double near_adc;
    int near_index;
    if(near_chan.X()!=-1){
      near_index = ReturnIndex(near_chan.X());
      near_adc = data->ADC(near_chan.X()-1);
      xtalk_adc[near_index]->Fill(near_adc);
      xtalk_rate = (near_adc - noise_mean[near_index]) / (cen_adc - noise_mean[trk_index]);
      if(xtalk_rate<0) xtalk_rate = 0;
      xtalk_esti[near_index]->Fill(xtalk_rate);
      xtalk_all->Fill(xtalk_rate);
    }
    if(near_chan.Y()!=-1){
      near_index = ReturnIndex(near_chan.Y());
      near_adc = data->ADC(near_chan.Y()-1);
      xtalk_adc[near_index]->Fill(near_adc);
      xtalk_rate = (near_adc - noise_mean[near_index]) / (cen_adc - noise_mean[trk_index]);
      if(xtalk_rate<0) xtalk_rate = 0;
      xtalk_esti[near_index]->Fill(xtalk_rate);
      xtalk_all->Fill(xtalk_rate);
    }  
    if(near_chan.Z()!=-1){
      near_index = ReturnIndex(near_chan.Z());
      near_adc = data->ADC(near_chan.Z()-1);
      xtalk_adc[near_index]->Fill(near_adc);
      xtalk_rate = (near_adc - noise_mean[near_index]) / (cen_adc - noise_mean[trk_index]);
      if(xtalk_rate<0) xtalk_rate = 0;
      xtalk_esti[near_index]->Fill(xtalk_rate);
      xtalk_all->Fill(xtalk_rate);
    }
    if(near_chan.T()!=-1){
      near_index = ReturnIndex(near_chan.T());
      near_adc = data->ADC(near_chan.T()-1);
      xtalk_adc[near_index]->Fill(near_adc);
      xtalk_rate = (near_adc - noise_mean[near_index]) / (cen_adc - noise_mean[trk_index]);
      if(xtalk_rate<0) xtalk_rate = 0;
      xtalk_esti[near_index]->Fill(xtalk_rate);
      xtalk_all->Fill(xtalk_rate);
    }

  }
  
  TText *pl_name = new TText();
  pl_name->SetTextSize(0.04);
  TText *pl_mean = new TText();
  pl_mean->SetTextSize(0.03);
  TText *pl_rms = new TText();
  pl_rms->SetTextSize(0.03);

  double st_left = 0.55;
  double st_top = 0.7;
  pl_mean->SetTextSize(0.04);
  pl_rms->SetTextSize(0.04);
  
  TCanvas *c1 = new TCanvas("noise","noise",1200,1200);
  c1->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c1->cd(i+1);
    noise_esti[i]->Draw("hist");
    name.Form("Overall mean = %f ADC",noise_mean[i]);
    pl_mean->DrawTextNDC(st_left,st_top,name);
    name.Form("Overall RMS = %f ADC",noise_rms[i]);
    pl_rms->DrawTextNDC(st_left,st_top-0.07,name);
  }
  c1->Update();

  TCanvas *c2 = new TCanvas("xtalk_adc","xtalk_adc",1200,1200);
  c2->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c2->cd(i+1);
    xtalk_adc[i]->Draw("hist");
  }
  c2->Update();
  
  TCanvas *c3 = new TCanvas("track_adc","track_adc",1200,1200);
  c3->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c3->cd(i+1);
    track_adc[i]->Draw("hist");
  }
  c3->Update();
  
  TCanvas *c4 = new TCanvas("xtalk_esti","xtalk_esti",1200,1200);
  c4->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c4->cd(i+1);
    xtalk_esti[i]->Draw("hist");
  }
  c4->Update();
  
  TCanvas *c5 = new TCanvas("xtalk_all","xtalk_all",700,600);
  c5->SetLeftMargin(0.15);
  c5->cd();
  xtalk_all->SetLineWidth(2);
  xtalk_all->SetLineColor(kBlue);
  //xtalk_all->SetMarkerColor(kBlue);
  //xtalk_all->SetMarkerStyle(20);
  //xtalk_all->Draw("P E0");
  xtalk_all->Draw("hist");

  name.Form("Crosstalk fraction:");
  pl_name->DrawTextNDC(0.5,0.68,name);
  name.Form("Mean = %f",xtalk_all->GetMean());
  pl_mean->DrawTextNDC(0.5,0.61,name);
  gPad->SetGridx();
  gPad->SetGridy();
  c5->Update();
  
  TString prefix = "../../../plots/scintillator_cube/";
  TString type = "3dprintcubes_newboard/";;
  TString suffix;

  suffix = prefix + type + "noise_esti.png";
  c1->SaveAs(suffix);
  
  suffix = prefix + type + "xtalk_adc.png";
  c2->SaveAs(suffix);

  suffix = prefix + type + "track_adc.png";
  c3->SaveAs(suffix);
  
  suffix = prefix + type + "xtalk_esti.png";
  c4->SaveAs(suffix);
  
  suffix = prefix + type + "xtalk_all.png";
  c5->SaveAs(suffix);

  suffix = prefix + type + "xtalk_all.pdf";
  c5->SaveAs(suffix);

  TString fout_name = "../../results/CrosstalkAnalysis_3DPrintCube_NewBoard.root";

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();
 
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c5->Write();
  
  for(int i = 0; i < fChanNum; i++){
    noise_esti[i]->Write();
    xtalk_adc[i]->Write();
    track_adc[i]->Write();
    xtalk_esti[i]->Write();
  }
 
  xtalk_all->Write();

  fout->Close();  

}

void DrawMPPCChanADC(){

  bool swap = false;
  bool newswap = false;

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();
  
  std::cout << "Total number of events: " << n_event << std::endl;
  
  TString name;
  TH1D *MPPC_ly[fChanNum];
  
  for(int i = 0; i < fChanNum; i++){
    name.Form("MPPC_ly_chan%i",fChanOrder[i]);
    MPPC_ly[i] = new TH1D(name,name,80,0,4500);
    MPPC_ly[i]->GetXaxis()->SetTitle("ADC");
    MPPC_ly[i]->GetYaxis()->SetTitle("Number of events / bin");
    MPPC_ly[i]->SetTitle(name);
    MPPC_ly[i]->SetLineWidth(2);
    MPPC_ly[i]->SetLineColor(kBlue);
  }

  for(int n = 0; n < n_event; n++){

    data->GetMppc(n,swap,newswap);
    
    for(int i = 0; i < fChanNum; i++) MPPC_ly[i]->Fill(data->ADC(fChanOrder[i]-1)); 

  }

  gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("MPPC2D","MPPC2D",1200,1200);
  c1->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c1->cd(i+1);
    MPPC_ly[i]->Draw();
    gPad->SetLogy();
  }
  c1->Update();

  TString prefix = "../../../plots/scintillator_cube/";
  TString type = "3dprintcubes_newboard/";;
  TString suffix;

  suffix = prefix + type + "MPPC2D.png";
  c1->SaveAs(suffix);

  TString fout_name = "../../results/MPPCLightYield_3DPrintCube_NewBoard.root";

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();
  
  for(int i = 0; i < fChanNum; i++) MPPC_ly[i]->Write();
 
  c1->Write();

  fout->Close();

}

// ---------------------------
// Auxiliary functions
// ---------------------------

TLorentzVector GetNearbyChannel(int chan){

  TLorentzVector near_chan;
  
  // Order: X = left, Y = top, Z = right, T = down
  if(chan==15) near_chan.SetXYZT(-1,-1,17,21);
  else if(chan==17) near_chan.SetXYZT(15,-1,19,23);
  else if(chan==19) near_chan.SetXYZT(17,-1,-1,25);
  else if(chan==21) near_chan.SetXYZT(-1,15,23,27);
  else if(chan==23) near_chan.SetXYZT(21,17,25,29);
  else if(chan==25) near_chan.SetXYZT(23,19,-1,31);
  else if(chan==27) near_chan.SetXYZT(-1,21,29,-1);
  else if(chan==29) near_chan.SetXYZT(27,23,31,-1);
  else if(chan==31) near_chan.SetXYZT(29,25,-1,-1);
  else near_chan.SetXYZT(-1,-1,-1,-1);

  return near_chan;

}

std::tuple<bool,int> CheckEventCosmicVertical(Mppc *data){

  // Step 1: Find the position of highest ADC cube in 3*3 matrix
  int top_index;
  
  double adc_temp = 0;
  for(int i = 0; i < fChanNum; i++){
    if(data->ADC(fChanOrder[i]-1)>=adc_temp){
      adc_temp = data->ADC(fChanOrder[i]-1);
      top_index = i;
    }
  }
  
  TVector3 top_pos = GetTopChanPos(fChanOrder[top_index]);
  
  // Step 2: Find the position of highest ADC cubes in trigger cube two layers

  // Top layer
  int midd_chan_xz, midd_chan_yz;
  
  adc_temp = 0;
  for(int i = 0; i < 3; i++){
    if(data->ADC(f13TopChanOrder[i]-1)>=adc_temp){
      adc_temp = data->ADC(f13TopChanOrder[i]-1);
      midd_chan_xz = f13TopChanOrder[i];
    }
  }
  
  adc_temp = 0;
  for(int i = 0; i < 3; i++){
    if(data->ADC(f24TopChanOrder[i]-1)>=adc_temp){
      adc_temp = data->ADC(f24TopChanOrder[i]-1);
      midd_chan_yz = f24TopChanOrder[i];
    }
  }

  TVector3 midd_pos = GetTrigChanPos(midd_chan_xz,midd_chan_yz);
  
  // Bottom layer
  int down_chan_xz, down_chan_yz;
  
  adc_temp = 0;
  for(int i = 0; i < 3; i++){
    if(data->ADC(f13DownChanOrder[i]-1)>=adc_temp){
      adc_temp = data->ADC(f13DownChanOrder[i]-1);
      down_chan_xz = f13DownChanOrder[i];
    }
  }
  
  adc_temp = 0;
  for(int i = 0; i < 3; i++){
    if(data->ADC(f24DownChanOrder[i]-1)>=adc_temp){
      adc_temp = data->ADC(f24DownChanOrder[i]-1);
      down_chan_yz = f24DownChanOrder[i];
    }
  }

  TVector3 down_pos = GetTrigChanPos(down_chan_xz,down_chan_yz);

  // Step 4: Check the track is vertical
  bool ver_tag;
  
  if(top_pos.X()==midd_pos.X() && midd_pos.X()==down_pos.X() &&
     top_pos.Y()==midd_pos.Y() && midd_pos.Y()==down_pos.Y()) ver_tag = true;
  else ver_tag = false;
  
  // Step 5: The trigger channel ADC should be larger than a threshold
  bool adc_tag;
  
  if(data->ADC(midd_chan_xz-1)>fADCCut &&
     data->ADC(midd_chan_yz-1)>fADCCut &&
     data->ADC(down_chan_xz-1)>fADCCut &&
     data->ADC(down_chan_yz-1)>fADCCut) adc_tag = true;
  else adc_tag = false;
  
  if(ver_tag==true && adc_tag==true) return std::make_tuple(true,top_index);
  else return std::make_tuple(false,-1);

}

TVector3 GetTopChanPos(int chan){

  TVector3 pos;
  
  if(chan==27) pos.SetXYZ(0,0,3);
  else if(chan==29) pos.SetXYZ(1,0,3);
  else if(chan==31) pos.SetXYZ(2,0,3);
  else if(chan==21) pos.SetXYZ(0,1,3);
  else if(chan==23) pos.SetXYZ(1,1,3);
  else if(chan==25) pos.SetXYZ(2,1,3);
  else if(chan==15) pos.SetXYZ(0,2,3);
  else if(chan==17) pos.SetXYZ(1,2,3);
  else if(chan==19) pos.SetXYZ(2,2,3);
  else pos.SetXYZ(-1,-1,-1);

  return pos;

}

TVector3 GetTrigChanPos(int x, int y){

  TVector3 pos;
  
  if(x==3 && y==11) pos.SetXYZ(0,0,0);
  else if(x==3 && y==13) pos.SetXYZ(0,1,0);
  else if(x==3 && y==7) pos.SetXYZ(0,2,0);
  else if(x==1 && y==11) pos.SetXYZ(1,0,0);
  else if(x==1 && y==13) pos.SetXYZ(1,1,0);
  else if(x==1 && y==7) pos.SetXYZ(1,2,0);
  else if(x==5 && y==11) pos.SetXYZ(2,0,0);
  else if(x==5 && y==13) pos.SetXYZ(2,1,0);
  else if(x==5 && y==7) pos.SetXYZ(2,2,0);
  else if(x==4 && y==12) pos.SetXYZ(0,0,2);
  else if(x==4 && y==14) pos.SetXYZ(0,1,2);
  else if(x==4 && y==8) pos.SetXYZ(0,2,2);
  else if(x==2 && y==12) pos.SetXYZ(1,0,2);
  else if(x==2 && y==14) pos.SetXYZ(1,1,2);
  else if(x==2 && y==8) pos.SetXYZ(1,2,2);
  else if(x==6 && y==12) pos.SetXYZ(2,0,2);
  else if(x==6 && y==14) pos.SetXYZ(2,1,2);
  else if(x==6 && y==8) pos.SetXYZ(2,2,2);
  else pos.SetXYZ(-1,-1,-1);

  return pos;

}

int ReturnIndex(int chan){

  int index;
  
  if(chan==15) index = 0;
  else if(chan==17) index = 1;
  else if(chan==19) index = 2;
  else if(chan==21) index = 3;
  else if(chan==23) index = 4;
  else if(chan==25) index = 5;
  else if(chan==27) index = 6;
  else if(chan==29) index = 7;
  else if(chan==31) index = 8;

  return index;

}
