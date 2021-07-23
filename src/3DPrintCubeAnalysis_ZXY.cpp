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
std::string fin_name = "../../inputs/3dprinted_cubes/3DPrintCubes_NewBoard_old.root";
//std::string fin_name = "../../inputs/3dprinted_cubes/3DMatrix_th280DAC_OR32_18June2021.root";

const int fChanNum = 9; // 3*3 matrix
int fChanOrder[fChanNum] = {15,17,19,21,23,25,27,29,31}; // 3*3 matrix channel number

double fChanGain[fChanNum] = {35.6556,36.9649,32.6565,30.7933,38.5202,36.4342,36.7935,36.7,37.0}; 
//double fChanGain[fChanNum] = {35.6556,36.9649,32.6565,30.7933,38.5202,36.4342,36.7935,35.4026,35.4026};

double f13TopChanOrder[3] = {8,14,12}; // Channel number for trigger cubes
double f13DownChanOrder[3] = {7,13,11};
double f24TopChanOrder[3] = {6,2,4};
double f24DownChanOrder[3] = {5,1,3};

const int fTriChanNum = 12;
int fTriChanOrder[fTriChanNum] = {8,14,12,7,13,11,6,2,4,5,1,3};

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
  TH1D *track_pe[fChanNum];
  
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
    xtalk_esti[i] = new TH1D(name,name,30,0,3);
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
    
    name.Form("channel%i_trackpe",fChanOrder[i]);
    track_pe[i] = new TH1D(name,name,50,0,120);
    track_pe[i]->GetXaxis()->SetTitle("Track channel light yield / p.e.");
    track_pe[i]->GetYaxis()->SetTitle("Number of events / bin");
    track_pe[i]->GetXaxis()->SetLabelSize(0.04);
    track_pe[i]->GetXaxis()->SetTitleSize(0.04);
    track_pe[i]->GetYaxis()->SetLabelSize(0.04);
    track_pe[i]->GetYaxis()->SetTitleSize(0.04);
    track_pe[i]->GetYaxis()->SetTitleOffset(1.4);
    track_pe[i]->SetTitle(name);
    track_pe[i]->SetLineWidth(2);
    track_pe[i]->SetLineColor(kBlue);    
  }
 
  TH1D *xtalk_all = new TH1D("xtalk_rate","xtalk_rate",30,0,3);
  xtalk_all->GetXaxis()->SetTitle("Crosstalk fraction / %");
  xtalk_all->GetYaxis()->SetTitle("Number of events / bin");
  xtalk_all->GetXaxis()->SetLabelSize(0.05);
  xtalk_all->GetXaxis()->SetTitleSize(0.05);
  xtalk_all->GetYaxis()->SetLabelSize(0.05);
  xtalk_all->GetYaxis()->SetTitleSize(0.05);
  xtalk_all->GetYaxis()->SetTitleOffset(1.4);
  xtalk_all->SetTitle("");
 
  TH1D *xtalk_allex = new TH1D("xtalk_rateex","xtalk_rateex",60,-3,3);
  xtalk_allex->GetXaxis()->SetTitle("Crosstalk fraction / %");
  xtalk_allex->GetYaxis()->SetTitle("Number of events / bin");
  xtalk_allex->GetXaxis()->SetLabelSize(0.05);
  xtalk_allex->GetXaxis()->SetTitleSize(0.05);
  xtalk_allex->GetYaxis()->SetLabelSize(0.05);
  xtalk_allex->GetYaxis()->SetTitleSize(0.05);
  xtalk_allex->GetYaxis()->SetTitleOffset(1.4);
  xtalk_allex->SetTitle("");

  TH1D *trkpe_all = new TH1D("trkpe_all","trkpe_all",50,0,120);
  trkpe_all->GetXaxis()->SetTitle("Track channel light yield / p.e.");
  trkpe_all->GetYaxis()->SetTitle("Number of events / bin");
  trkpe_all->GetXaxis()->SetLabelSize(0.05);
  trkpe_all->GetXaxis()->SetTitleSize(0.05);
  trkpe_all->GetYaxis()->SetLabelSize(0.05);
  trkpe_all->GetYaxis()->SetTitleSize(0.05);
  trkpe_all->GetYaxis()->SetTitleOffset(1.4);
  trkpe_all->SetTitle("");
  
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
    double cen_gain = fChanGain[trk_index];
    if(cen_adc<400) continue;
    track_adc[trk_index]->Fill(cen_adc);
    double cen_temp = (cen_adc - noise_mean[trk_index]) / cen_gain;
    track_pe[trk_index]->Fill(cen_temp);
    trkpe_all->Fill(cen_temp);
    
    TLorentzVector near_chan = GetNearbyChannel(fChanOrder[trk_index]);
    
    double xtalk_rate;
    double near_adc;
    double near_gain;
    int near_index;
    if(near_chan.X()!=-1){
      near_index = ReturnIndex(near_chan.X());
      near_adc = data->ADC(near_chan.X()-1);
      near_gain = fChanGain[near_index];
      xtalk_adc[near_index]->Fill(near_adc);
      xtalk_rate = ((near_adc - noise_mean[near_index]) / near_gain) / ((cen_adc - noise_mean[trk_index]) / cen_gain);
      xtalk_allex->Fill(xtalk_rate);
      if(xtalk_rate<0) xtalk_rate = 0;
      xtalk_esti[near_index]->Fill(xtalk_rate);
      xtalk_all->Fill(xtalk_rate);
    }
    if(near_chan.Y()!=-1){
      near_index = ReturnIndex(near_chan.Y());
      near_adc = data->ADC(near_chan.Y()-1);
      near_gain = fChanGain[near_index];
      xtalk_adc[near_index]->Fill(near_adc);
      xtalk_rate = ((near_adc - noise_mean[near_index]) / near_gain) / ((cen_adc - noise_mean[trk_index]) / cen_gain);
      xtalk_allex->Fill(xtalk_rate);
      if(xtalk_rate<0) xtalk_rate = 0;
      xtalk_esti[near_index]->Fill(xtalk_rate);
      xtalk_all->Fill(xtalk_rate);
    }  
    if(near_chan.Z()!=-1){
      near_index = ReturnIndex(near_chan.Z());
      near_adc = data->ADC(near_chan.Z()-1);
      near_gain = fChanGain[near_index];
      xtalk_adc[near_index]->Fill(near_adc);
      xtalk_rate = ((near_adc - noise_mean[near_index]) / near_gain) / ((cen_adc - noise_mean[trk_index]) / cen_gain);
      xtalk_allex->Fill(xtalk_rate);
      if(xtalk_rate<0) xtalk_rate = 0;
      xtalk_esti[near_index]->Fill(xtalk_rate);
      xtalk_all->Fill(xtalk_rate);
    }
    if(near_chan.T()!=-1){
      near_index = ReturnIndex(near_chan.T());
      near_adc = data->ADC(near_chan.T()-1);
      near_gain = fChanGain[near_index];
      xtalk_adc[near_index]->Fill(near_adc);
      xtalk_rate = ((near_adc - noise_mean[near_index]) / near_gain) / ((cen_adc - noise_mean[trk_index]) / cen_gain);
      xtalk_allex->Fill(xtalk_rate);
      if(xtalk_rate<0) xtalk_rate = 0;
      xtalk_esti[near_index]->Fill(xtalk_rate);
      xtalk_all->Fill(xtalk_rate);
    }

  }
 
  // Compute mean light yield per channel
  double ly_mean[fChanNum];
  double ly_overall = trkpe_all->GetMean();
  for(int i = 0; i < fChanNum; i++) ly_mean[i] = track_pe[i]->GetMean();

  TH1D *chanly_ave = new TH1D("chanly_ave","chanly_ave",fChanNum,0,9);
  chanly_ave->GetXaxis()->SetTitle("Channel");
  chanly_ave->GetYaxis()->SetTitle("Average light yield / p.e.");
  chanly_ave->GetXaxis()->SetLabelSize(0.05);
  chanly_ave->GetXaxis()->SetTitleSize(0.05);
  chanly_ave->GetYaxis()->SetLabelSize(0.05);
  chanly_ave->GetYaxis()->SetTitleSize(0.05);
  chanly_ave->GetYaxis()->SetTitleOffset(1.4);
  chanly_ave->SetTitle("");

  TGraph *ly_general = new TGraph(fChanNum);
  
  TString label;
  for(int i = 0; i < fChanNum; i++){
    label.Form("%i",fChanOrder[i]);
    chanly_ave->GetXaxis()->SetBinLabel(i+1,label);
    chanly_ave->SetBinContent(i+1,ly_mean[i]);
    ly_general->SetPoint(i,i,ly_overall);
  }
  ly_general->SetPoint(fChanNum,9,ly_overall);
  
  chanly_ave->SetLineColor(kBlue);
  chanly_ave->SetLineWidth(2);
  chanly_ave->GetYaxis()->SetRangeUser(0,80);
  ly_general->SetLineColor(kRed);
  ly_general->SetLineWidth(2);
  ly_general->SetLineStyle(2);
  
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
  
  gStyle->SetOptStat(0);
  
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
    gPad->SetLogy();
  }
  c4->Update();
  
  TCanvas *c5 = new TCanvas("xtalk_all","xtalk_all",700,600);
  c5->SetLeftMargin(0.15);
  c5->cd();
  xtalk_all->SetLineWidth(2);
  xtalk_all->SetLineColor(kRed);
  xtalk_all->Draw("hist");

  name.Form("Crosstalk fraction:");
  pl_name->DrawTextNDC(0.5,0.68,name);
  name.Form("Mean = %f",xtalk_all->GetMean());
  pl_mean->DrawTextNDC(0.5,0.61,name);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetTicks(1,1);
  c5->Update();

  TCanvas *c6 = new TCanvas("xtalk_allex","xtalk_allex",700,600);
  c6->SetLeftMargin(0.15);
  c6->cd();
  xtalk_allex->SetLineWidth(2);
  xtalk_allex->SetLineColor(kRed);
  xtalk_allex->Draw("hist");

  name.Form("Crosstalk fraction:");
  pl_name->DrawTextNDC(0.5,0.68,name);
  name.Form("Mean = %f",xtalk_allex->GetMean());
  pl_mean->DrawTextNDC(0.5,0.61,name);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  gPad->SetTicks(1,1);
  c6->Update();
 
  TCanvas *c7 = new TCanvas("track_pe","track_pe",1200,1200);
  c7->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c7->cd(i+1);
    track_pe[i]->Draw("hist");
    name.Form("Overall mean = %f p.e.",ly_mean[i]);
    pl_mean->DrawTextNDC(st_left,st_top,name);
  }
  c7->Update();
  
  TCanvas *c8 = new TCanvas("trkpe_all","trkpe_all",700,600);
  c8->SetLeftMargin(0.15);
  c8->cd();
  trkpe_all->SetLineWidth(2);
  trkpe_all->SetLineColor(kRed);
  trkpe_all->Draw("hist");

  //name.Form("Crosstalk fraction:");
  //pl_name->DrawTextNDC(0.5,0.68,name);
  name.Form("Overall mean = %f p.e.",ly_overall);
  pl_mean->DrawTextNDC(0.5,0.61,name);
  gPad->SetGridx();
  gPad->SetGridy();
  //gPad->SetLogy();
  gPad->SetTicks(1,1);
  c8->Update();  
  
  TCanvas *c9 = new TCanvas("chanly_comp","chanly_comp",700,600);
  c9->SetLeftMargin(0.15);
  c9->cd();
  chanly_ave->Draw("hist");
  ly_general->Draw("l");
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetTicks(1,1);
  c9->Update();
  
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

  suffix = prefix + type + "xtalk_allex.png";
  c6->SaveAs(suffix);

  suffix = prefix + type + "xtalk_allex.pdf";
  c6->SaveAs(suffix);

  suffix = prefix + type + "track_pe.png";
  c7->SaveAs(suffix);

  suffix = prefix + type + "trkpe_all.png";
  c8->SaveAs(suffix);

  suffix = prefix + type + "trkpe_all.pdf";
  c8->SaveAs(suffix);

  suffix = prefix + type + "chanly_comp.png";
  c9->SaveAs(suffix);

  suffix = prefix + type + "chanly_comp.pdf";
  c9->SaveAs(suffix);

  TString fout_name = "../../results/CrosstalkAnalysis_3DPrintCube_NewBoard.root";

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();
 
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c5->Write();
  c6->Write();
  c7->Write();
  c8->Write();
  c9->Write();
  
  for(int i = 0; i < fChanNum; i++){
    noise_esti[i]->Write();
    xtalk_adc[i]->Write();
    track_adc[i]->Write();
    xtalk_esti[i]->Write();
    track_pe[i]->Write();
  }
 
  xtalk_all->Write();
  xtalk_allex->Write();
  trkpe_all->Write();
  chanly_ave->Write();
  ly_general->Write();

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
  TH1D *MPPCtri_ly[fTriChanNum];
  
  for(int i = 0; i < fChanNum; i++){
    name.Form("MPPC_ly_chan%i",fChanOrder[i]);
    MPPC_ly[i] = new TH1D(name,name,80,0,4500);
    MPPC_ly[i]->GetXaxis()->SetTitle("ADC");
    MPPC_ly[i]->GetYaxis()->SetTitle("Number of events / bin");
    MPPC_ly[i]->SetTitle(name);
    MPPC_ly[i]->SetLineWidth(2);
    MPPC_ly[i]->SetLineColor(kBlue);
  }

  for(int i = 0; i < fTriChanNum; i++){
    name.Form("MPPC_ly_chan%i",fTriChanOrder[i]);
    MPPCtri_ly[i] = new TH1D(name,name,80,0,4500);
    MPPCtri_ly[i]->GetXaxis()->SetTitle("ADC");
    MPPCtri_ly[i]->GetYaxis()->SetTitle("Number of events / bin");
    MPPCtri_ly[i]->SetTitle(name);
    MPPCtri_ly[i]->SetLineWidth(2);
    MPPCtri_ly[i]->SetLineColor(kBlue);
  }

  for(int n = 0; n < n_event; n++){

    data->GetMppc(n,swap,newswap);
    
    for(int i = 0; i < fChanNum; i++) MPPC_ly[i]->Fill(data->ADC(fChanOrder[i]-1)); 
    
    for(int i = 0; i < fTriChanNum; i++) MPPCtri_ly[i]->Fill(data->ADC(fTriChanOrder[i]-1));

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

  TCanvas *c2 = new TCanvas("MPPC2D_13","MPPC2D_13",1200,800);
  c2->Divide(3,2);
  for(int i = 0; i < 6; i++){
    c2->cd(i+1);
    MPPCtri_ly[i]->Draw();
    gPad->SetLogy();
  }
  c2->Update();

  TCanvas *c3 = new TCanvas("MPPC2D_24","MPPC2D_24",1200,800);
  c3->Divide(3,2);
  for(int i = 0; i < 6; i++){
    c3->cd(i+1);
    MPPCtri_ly[i+6]->Draw();
    gPad->SetLogy();
  }
  c3->Update();

  TString prefix = "../../../plots/scintillator_cube/";
  TString type = "3dprintcubes_newboard/";;
  TString suffix;

  suffix = prefix + type + "MPPC2D.png";
  c1->SaveAs(suffix);

  suffix = prefix + type + "MPPC2D_13.png";
  c2->SaveAs(suffix);
  
  suffix = prefix + type + "MPPC2D_24.png";
  c3->SaveAs(suffix);

  TString fout_name = "../../results/MPPCLightYield_3DPrintCube_NewBoard.root";

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();
  
  for(int i = 0; i < fChanNum; i++) MPPC_ly[i]->Write();
  
  for(int i = 0; i < fTriChanNum; i++) MPPCtri_ly[i]->Write();
 
  c1->Write();
  c2->Write();
  c3->Write();

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
  
  if(/*top_pos.X()==midd_pos.X() &&*/ midd_pos.X()==down_pos.X() &&
     /*top_pos.Y()==midd_pos.Y() &&*/ midd_pos.Y()==down_pos.Y()) ver_tag = true;
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
 
  if(x==7 && y==5) pos.SetXYZ(0,0,0);
  else if(x==7 && y==1) pos.SetXYZ(0,1,0);
  else if(x==7 && y==3) pos.SetXYZ(0,2,0);
  else if(x==13 && y==5) pos.SetXYZ(1,0,0);
  else if(x==13 && y==1) pos.SetXYZ(1,1,0);
  else if(x==13 && y==3) pos.SetXYZ(1,2,0);
  else if(x==11 && y==5) pos.SetXYZ(2,0,0);
  else if(x==11 && y==1) pos.SetXYZ(2,1,0);
  else if(x==11 && y==3) pos.SetXYZ(2,2,0);
  else if(x==8 && y==6) pos.SetXYZ(0,0,2);
  else if(x==8 && y==2) pos.SetXYZ(0,1,2);
  else if(x==8 && y==4) pos.SetXYZ(0,2,2);
  else if(x==14 && y==6) pos.SetXYZ(1,0,2);
  else if(x==14 && y==2) pos.SetXYZ(1,1,2);
  else if(x==14 && y==4) pos.SetXYZ(1,2,2);
  else if(x==12 && y==6) pos.SetXYZ(2,0,2);
  else if(x==12 && y==2) pos.SetXYZ(2,1,2);
  else if(x==12 && y==4) pos.SetXYZ(2,2,2); 
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
