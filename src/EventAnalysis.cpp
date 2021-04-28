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
bool CheckIndexMatch(int,int,int,int);
bool EventCosmicCut(std::vector<double>,double);

// Draw event 3D display (after some cuts)
// Notation: file_option = input data file want to be used
// Notation: seed = random number generator
// Notation: ADC_cut = a cut applied on ADC value
void DrawEvent3D(int file_option = 1, int seed = 0, double ADC_cut = 500){

  std::string fin_name;

  if(file_option==1){
    fin_name = "../SpecialRun_GluedCubesNoTeflon_14April2021.root";
  }
  else if(file_option==2){
    fin_name = "../SpecialRun_GluedCubesNoTeflon_15April2021.root";
  }
  else if(file_option==3){
    fin_name = "../SpecialRun_GluedCubesNoTeflon_16April2021.root";
  }
  else if(file_option==4){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_19April2021.root";
  }
  else if(file_option==5){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_19April2021_2nd.root";
  }
  else if(file_option==6){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_20April2021.root";
  }
  else if(file_option==7){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_20April2021_2nd.root";
  }
  else if(file_option==8){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_21April2021.root";
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

      // In order to draw 3D event, at each layer only consider the highest ADC on each plane

      // Top layer (XZ plane)
      ADC_max = 0;
      if(data->ADC(13)>=ADC_max){
        ADC_max = data->ADC(13);
        index_xz = 14;
        ADC_xz = data->ADC(13);
      }
      if(data->ADC(3)>=ADC_max){
        ADC_max = data->ADC(3);
        index_xz = 4;
        ADC_xz = data->ADC(3);
      }
      if(data->ADC(12)>=ADC_max){
        ADC_max = data->ADC(12);
        index_xz = 13;
        ADC_xz = data->ADC(12);
      }
      // Top layer (YZ plane)
      ADC_max = 0;
      if(data->ADC(8)>=ADC_max){ 
        ADC_max = data->ADC(8);
        index_yz = 9;
        ADC_yz = data->ADC(8);
      }
      if(data->ADC(17)>=ADC_max){ 
        ADC_max = data->ADC(17);
        index_yz = 18;
        ADC_yz = data->ADC(17);
      }
      if(data->ADC(7)>=ADC_max){ 
        ADC_max = data->ADC(7);
        index_yz = 8;
        ADC_yz = data->ADC(7);
      }
      
      cube_pos = GetMatchCubePos(index_xz,index_yz);
      Event_show->Fill(cube_pos.X(),cube_pos.Y(),cube_pos.Z(),ADC_xz+ADC_yz); 

      // Middle layer (XZ plane)
      ADC_max = 0;
      if(data->ADC(2)>=ADC_max){
        ADC_max = data->ADC(2);
        index_xz = 3;
        ADC_xz = data->ADC(2);
      }
      if(data->ADC(11)>=ADC_max){
        ADC_max = data->ADC(11);
        index_xz = 12;
        ADC_xz = data->ADC(11);
      }
      if(data->ADC(1)>=ADC_max){
        ADC_max = data->ADC(1);
        index_xz = 2;
        ADC_xz = data->ADC(1);
      }
      // Middle layer (YZ plane)
      ADC_max = 0;
      if(data->ADC(16)>=ADC_max){
        ADC_max = data->ADC(16);
        index_yz = 17;
        ADC_yz = data->ADC(16);
      }
      if(data->ADC(6)>=ADC_max){
        ADC_max = data->ADC(6);
        index_yz = 7;
        ADC_yz = data->ADC(6);
      }
      if(data->ADC(15)>=ADC_max){
        ADC_max = data->ADC(15);
        index_yz = 16;
        ADC_yz = data->ADC(15);
      }

      cube_pos = GetMatchCubePos(index_xz,index_yz);
      Event_show->Fill(cube_pos.X(),cube_pos.Y(),cube_pos.Z(),ADC_xz+ADC_yz);

      // Bottom layer (XZ plane)
      ADC_max = 0;
      if(data->ADC(10)>=ADC_max){
        ADC_max = data->ADC(10);
        index_xz = 11;
        ADC_xz = data->ADC(10);
      }
      if(data->ADC(0)>=ADC_max){
        ADC_max = data->ADC(0);
        index_xz = 1;
        ADC_xz = data->ADC(0);
      }
      if(data->ADC(9)>=ADC_max){
        ADC_max = data->ADC(9);
        index_xz = 10;
        ADC_xz = data->ADC(9);
      }
      // Bottom layer (YZ plane)
      ADC_max = 0;
      if(data->ADC(5)>=ADC_max){
        ADC_max = data->ADC(5);
        index_yz = 6;
        ADC_yz = data->ADC(5);
      }
      if(data->ADC(14)>=ADC_max){
        ADC_max = data->ADC(14);
        index_yz = 15;
        ADC_yz = data->ADC(14);
      }
      if(data->ADC(4)>=ADC_max){
        ADC_max = data->ADC(4);
        index_yz = 5;
        ADC_yz = data->ADC(4);
      }

      cube_pos = GetMatchCubePos(index_xz,index_yz);
      Event_show->Fill(cube_pos.X(),cube_pos.Y(),cube_pos.Z(),ADC_xz+ADC_yz);

      // Also draw 2D MPPC on each plane
      for(int i = 0; i < 18; i++){
        cube_pos = GetMPPC2DPos(i+1);
        MPPC2D_xz->Fill(cube_pos.X(),cube_pos.Z(),data->ADC(i));
        MPPC2D_yz->Fill(cube_pos.Y(),cube_pos.Z(),data->ADC(i)); 
      }
 
      out_tag = true; 
   
    }

  } 

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
  c4->Update();

  // Save the plot into output file
  TFile *fout = new TFile("./EventDisplay_AfterADCCut.root","recreate");
  fout->cd();

  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();

  Event_show->Write();
  MPPC2D_xz->Write();
  MPPC2D_yz->Write();

  fout->Close();

}

// Plot the ADC distribution for each channel
void DrawMPPCLightYield(int file_option = 1, double ADC_cut = 500){

  std::string fin_name;

  if(file_option==1){
    fin_name = "../SpecialRun_GluedCubesNoTeflon_14April2021.root";
  }
  else if(file_option==2){
    fin_name = "../SpecialRun_GluedCubesNoTeflon_15April2021.root";
  }
  else if(file_option==3){
    fin_name = "../SpecialRun_GluedCubesNoTeflon_16April2021.root";
  }
  else if(file_option==4){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_19April2021.root";
  }
  else if(file_option==5){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_19April2021_2nd.root";
  }
  else if(file_option==6){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_20April2021.root";
  }
  else if(file_option==7){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_20April2021_2nd.root";
  }
  else if(file_option==8){
    fin_name = "../SpecialRun_GluedCubesWithTeflon_21April2021.root";
  }

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();

  // Light yield distributions for each channel 0 - 17
  TString name;
  TH1D *MPPC_ly[18];
  
  for(int i = 0; i < 18; i++){
    name.Form("MPPC_ly_chan%i",i+1);
    MPPC_ly[i] = new TH1D(name,name,100,0,4000);
    MPPC_ly[i]->GetXaxis()->SetTitle("ADC");
    MPPC_ly[i]->GetYaxis()->SetTitle("Number of events / bin");
    MPPC_ly[i]->SetTitle(name);
    MPPC_ly[i]->SetLineWidth(2);
    MPPC_ly[i]->SetLineColor(kBlue);
  }

  // Loop over all events
  for(int n = 0; n < n_event; n++){

    data->GetMppc(n);

    for(int i = 0; i < 18; i++){
      if((data->ADC(i))>ADC_cut) MPPC_ly[i]->Fill(data->ADC(i)); 
    }

  }

  gStyle->SetOptStat(0);

  // MPPC channel face 1 + 3
  TCanvas *c1 = new TCanvas("MPPC2D_xz","MPPC2D_xz",1200,1200);
  c1->Divide(3,3);
  c1->cd(1);
  MPPC_ly[13]->Draw();
  gPad->SetLogy();
  c1->cd(2);
  MPPC_ly[4]->Draw();
  gPad->SetLogy();
  c1->cd(3);
  MPPC_ly[12]->Draw();
  gPad->SetLogy();
  c1->cd(4);
  MPPC_ly[2]->Draw();
  gPad->SetLogy(); 
  c1->cd(5);
  MPPC_ly[11]->Draw();
  gPad->SetLogy();
  c1->cd(6);
  MPPC_ly[1]->Draw();
  gPad->SetLogy();
  c1->cd(7);
  MPPC_ly[10]->Draw();
  gPad->SetLogy();
  c1->cd(8);
  MPPC_ly[0]->Draw();
  gPad->SetLogy();
  c1->cd(9);
  MPPC_ly[9]->Draw();
  gPad->SetLogy();
  c1->Update();

  // MPPC channel face 2 + 4
  TCanvas *c2 = new TCanvas("MPPC2D_yz","MPPC2D_yz",1200,1200);
  c2->Divide(3,3);
  c2->cd(1);
  MPPC_ly[8]->Draw();
  gPad->SetLogy();
  c2->cd(2);
  MPPC_ly[17]->Draw();
  gPad->SetLogy();
  c2->cd(3);
  MPPC_ly[7]->Draw(); 
  gPad->SetLogy();
  c2->cd(4);
  MPPC_ly[16]->Draw();
  gPad->SetLogy();
  c2->cd(5);
  MPPC_ly[6]->Draw();
  gPad->SetLogy();
  c2->cd(6);
  MPPC_ly[15]->Draw();
  gPad->SetLogy();
  c2->cd(7);
  MPPC_ly[5]->Draw();
  gPad->SetLogy();
  c2->cd(8);
  MPPC_ly[14]->Draw();
  gPad->SetLogy();
  c2->cd(9);
  MPPC_ly[4]->Draw();
  gPad->SetLogy();
  c2->Update();

  // Save the plots into output file
  TFile *fout = new TFile("./MPPC_ChannelLightYield.root","recreate");
  fout->cd();
  
  for(int i = 0; i < 18; i++) MPPC_ly[i]->Write();
 
  c1->Write();
  c2->Write();

  fout->Close();

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
