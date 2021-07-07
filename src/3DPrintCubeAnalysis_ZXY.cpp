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

// Global parameters
std::string fin_name = "../../inputs/3dprinted_cubes/3DPrintCubes_NewBoard.root";

const int fChanNum = 9; // 3*3 matrix
int fChanOrder[fChanNum] = {15,17,19,21,23,25,27,29,31};

// ---------------------------
// Main functions
// ---------------------------

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
