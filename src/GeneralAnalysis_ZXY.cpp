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

// ----------------------------------------

// Auxiliary function prototypes
TVector3 GetMPPC2DPos(int);
TVector3 GetMatchCubePos(int,int);
std::vector<TLorentzVector> Get3DMatchedCubes(std::vector<double>,double);
std::tuple<std::vector<double>,std::vector<double>,TGraph2D*> Get3DTrackFit(std::vector<TLorentzVector>);
bool CheckIndexMatch(int,int,int,int);
bool EventCosmicCut(std::vector<double>,double);
TVector3 GetCubeChannel(TVector3);
TLorentzVector GetNearbyChannel(int);
int ReturnIndex(int);
bool CheckTrackVertical(std::vector<TLorentzVector>);
template <class T> void SetHistoSettings(T*);

// ----------------------------------------

// Global parameters
// Cube size
const double fCubeSize = 10; // In mm
// MPPC channel numbers
const int fChanNum = 18; // MPPC channels
const int fChanTotNum = 32; // Total MPPC channels
// The ADC cut, currently for glued cubes (> 500 units) and for SFGD cubes (> 1000 units)
const double fADCCut = 1000; // units
// Maximum ADC cut, currently chosen to be 4000 ADC
const double fADCUppCut = 4000; // units

// ----------------------------------------

// Number of input files
const int fFileNum = 4;

// File names array
std::string fFileName[fFileNum] = {"glued_cubes/GluedCubes_OldBoard", // May 18, 19 and June 4
"sfgd_cubes/SFGDCubes_OldBoard", // May 4 - 7
"sfgd_cubes/SFGDCubes_NewBoard", // June 5, 7 - 9
"glued_cubes/GluedCubes_NewBoard" // June 10, 11
};

// Output file names array
std::string fOutFileName[fFileNum] = {"GluedCubes_OldBoard",
"SFGDCubes_OldBoard",
"SFGDCubes_NewBoard",
"GluedCubes_NewBoard"};

// Path names array
TString fPathName[fFileNum] = {"gluedcubes_oldboard",
"sfgdcubes_oldboard",
"sfgdcubes_newboard",
"gluedcubes_newboard"};

// Swap index
// Only SFGD cube data (old board) is swapped
bool fSwap[fFileNum] = {false,true,false,false};
// Only SFGD cube data (new board) is swapped
bool fNewSwap[fFileNum] = {false,false,true,false};

// Channel mapping choice, for new board, the channel is different from old board
// 1 = old board, 2 = new board
int fChanMapChoice;

// Drawing order array
// First index = 0 (old board), 1 (new board)
int fFace13[2][9] = {{13,4,14,3,12,2,10,1,11},{4,2,6,17,19,15,3,1,5}};
int fFace24[2][9] = {{9,18,8,16,7,17,6,15,5},{12,14,8,23,21,25,11,13,7}};

// Channel order
// First index = 0 (old board), 1 (new board)
int fChanOrder[2][18] = {

// Old mapping
{13,4,14,3,12,2,10,1,11,9,18,8,16,7,17,6,15,5},

// New mapping
{4,2,6,17,19,15,3,1,5,12,14,8,23,21,25,11,13,7}

};

// MPPC Gain
// First index = 0 (glued cube), 1 (SFGD cube)
double fMPPCGain[2][18] = {

// Glued cubes (old board)
//{42.5371,39.4829,42.5123,41.0193,42.2430,43.0222,41.4896,40.4243,40.5144, // Face 1 + 3
//40.0000,41.6886,40.7666,36.7137,38.7410,42.7172,38.7927,37.6620,38.6067}, // Face 2 + 4

// Glued cubes (new board)
{40.3914,39.2136,48.9645,38.1267,36.4901,37.1654,40.0,37.7858,37.675, // Face 1 + 3
42.5401,38.2835,41.5478,54.1939,34.7,36.1043,34.8194,42.5503,37.0925}, // Face 2 + 4

// SFGD cubes (old board)
//{34.4433,39.7124,35.4863,39.2688,37.4597,38.0388,36.9028,40.3587,41.8002, // Face 1 + 3
//43.3295,34.2956,35.6027,36.7148,36.4768,37.3488,36.7135,41.0380,36.7928} // Face 2 + 4

// SFGD cubes (new board)
{39.533,37.871,38.012,39.381,35.342,37.211,41.4896,37.439,38.46, // Face 1 + 3
37.56,36.945,40.4,37.53,38.606,43.763,43.59,37.943,35.315} // Face 2 + 4

};

// ----------------------------------------

// General cube crosstalk analysis
void CrosstalkAnalysis(int file_option = 1, double ADC_cut = fADCCut){

  std::string fin_name;
  bool swap, newswap;

  // Set input file name
  fin_name = "../../inputs/" + fFileName[file_option-1] + ".root";
  swap = fSwap[file_option-1];
  newswap = fNewSwap[file_option-1];

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();
  
  std::cout << "Total number of events: " << n_event << std::endl;

  // Channel mapping choice
  if(file_option==1 || file_option==2) fChanMapChoice = 1;
  else if(file_option==3 || file_option==4) fChanMapChoice = 2;
  
  int chan_order[fChanNum];
  for(int i = 0; i < fChanNum; i++) chan_order[i] = fChanOrder[fChanMapChoice-1][i];

  // Gain array
  double ADCgain[fChanNum];
  for(int i = 0; i < fChanNum; i++){
    if(file_option==1 || file_option==4) ADCgain[i] = fMPPCGain[0][i]; // Glued cubes
    else if(file_option==2 || file_option==3) ADCgain[i] = fMPPCGain[1][i]; // SFGD cubes
  }

  // Create some variables
  std::vector<double> ADC_temp;
  bool pass_tag;
  std::vector<TLorentzVector> cube_array;
  std::vector<double> track_info;
  std::vector<double> chan_path;
  TGraph2D *trk_graph;
  TVector3 cubepos_cen;
  TVector3 cubepos_bei[4];
  TVector3 cubechan_cen;
  TVector3 cubechan_bei[4];
  int nx, ny;
  double cubely_bei[4];
  double cubenoise_bei[4];
  double cubegain_bei[4];
  double xtalk_frac;
  double noise_frac;
  double noise_temp;
  TLorentzVector channel_near;
  int index_temp;

  // Create some histograms
  // Overall crosstalk rate (with ADC)
  TH1D *xtalk_rate = new TH1D("xtalk_rate","xtalk_rate",30,0,15);
  xtalk_rate->GetXaxis()->SetTitle("Crosstalk fraction / %");
  xtalk_rate->GetYaxis()->SetTitle("Number of events / bin");
  xtalk_rate->GetXaxis()->SetLabelSize(0.05);
  xtalk_rate->GetXaxis()->SetTitleSize(0.05);
  xtalk_rate->GetYaxis()->SetLabelSize(0.05);
  xtalk_rate->GetYaxis()->SetTitleSize(0.05);
  xtalk_rate->GetYaxis()->SetTitleOffset(1.4);
  xtalk_rate->SetTitle("");

  TH1D *comp_xtalkfrac = (TH1D*)xtalk_rate->Clone("comp_xtalkfrac");
  TH1D *comp_noisefrac = (TH1D*)xtalk_rate->Clone("comp_noisefrac");

  // Overall crosstalk rate (with P.E. from gain)
  TH1D *xtalk_rate_pe = (TH1D*)xtalk_rate->Clone("xtalk_rate_pe");
  TH1D *comp_xtalkfrac_pe = (TH1D*)xtalk_rate->Clone("comp_xtalkfrac_pe");
  TH1D *comp_noisefrac_pe = (TH1D*)xtalk_rate->Clone("comp_noisefrac_pe");

  // Channel ADC
  TH1D *trkcube_ly = new TH1D("trkcube_ly","trkcube_ly",90,0,4500);
  trkcube_ly->GetXaxis()->SetTitle("MPPC channel value / ADC");
  trkcube_ly->GetYaxis()->SetTitle("Probability density / bin");
  trkcube_ly->GetXaxis()->SetLabelSize(0.04);
  trkcube_ly->GetXaxis()->SetTitleSize(0.04);
  trkcube_ly->GetYaxis()->SetLabelSize(0.04);
  trkcube_ly->GetYaxis()->SetTitleSize(0.04);
  trkcube_ly->GetYaxis()->SetTitleOffset(1.4);
  trkcube_ly->SetTitle("");

  TH1D *beicube_ly = (TH1D*)trkcube_ly->Clone("beicube_ly");
  TH1D *noise_ly = (TH1D*)trkcube_ly->Clone("noise_ly");

  /*TH1D *noise_esti = new TH1D("noise_esti","noise_esti",60,0,400);
  noise_esti->GetXaxis()->SetTitle("Estimated noise level / ADC");
  noise_esti->GetYaxis()->SetTitle("Number of events / bin");
  noise_esti->GetXaxis()->SetLabelSize(0.04);
  noise_esti->GetXaxis()->SetTitleSize(0.04);
  noise_esti->GetYaxis()->SetLabelSize(0.04);
  noise_esti->GetYaxis()->SetTitleSize(0.04);
  noise_esti->GetYaxis()->SetTitleOffset(1.4);
  noise_esti->SetTitle("");*/

  // Noise level per channel
  TH1D *noise_esti[fChanNum];
  // Crosstalk level per channel
  TH1D *xtalk_esti[fChanNum];
  // Crosstalk ADC per channel
  TH1D *xtalk_adc[fChanNum];
  // Track ADC per channel
  TH1D *track_adc[fChanNum];

  TString name;
  for(int i = 0; i < fChanNum; i++){
    name.Form("channel%i_noise",chan_order[i]);
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

    name.Form("channel%i_xtalk",chan_order[i]);
    xtalk_esti[i] = new TH1D(name,name,60,0,30);
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

    name.Form("channel%i_xtalkADC",chan_order[i]);
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

    name.Form("channel%i_trackADC",chan_order[i]);
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

  // 2D plot, xtalk rate vs. track channel p.e.
  TH2D *test_2d = new TH2D("test_2d","test_2d",60,0,120,40,0,12);
  test_2d->GetXaxis()->SetTitle("Track channel / p.e.");
  test_2d->GetYaxis()->SetTitle("Crosstalk channel / p.e.");
  test_2d->SetTitle("");

  // Estimate the noise level in order to subtract it
  std::cout << "Start to estimate noise level" << std::endl;
  // First loop over all events
  for(int n = 0; n < n_event; n++){

    data->GetMppc(n,swap,newswap);

    ADC_temp.clear();
    for(int i = 0; i < fChanNum; i++) ADC_temp.push_back(data->ADC(chan_order[i]-1));

    // Method 1: use an overall noise level
    /*pass_tag = EventCosmicCut(ADC_temp,ADC_cut);
    if(pass_tag==false) continue;

    cube_array = Get3DMatchedCubes(ADC_temp,ADC_cut);
 
    // Only select the vertical tracks
    std::tie(track_info,chan_path,trk_graph) = Get3DTrackFit(cube_array);
    if(TMath::Abs(track_info[1])>1e-4) continue;

    // In order to estimate the noise level consider the cube furtherest away from the main cube
    cubepos_cen.SetXYZ(cube_array[1].X(),cube_array[1].Y(),cube_array[1].Z());
    cubepos_bei[0].SetXYZ(cube_array[1].X()-2,cube_array[1].Y(),cube_array[1].Z());
    cubepos_bei[1].SetXYZ(cube_array[1].X()+2,cube_array[1].Y(),cube_array[1].Z());
    cubepos_bei[2].SetXYZ(cube_array[1].X(),cube_array[1].Y()+2,cube_array[1].Z());
    cubepos_bei[3].SetXYZ(cube_array[1].X(),cube_array[1].Y()-2,cube_array[1].Z());
 
    cubechan_cen = GetCubeChannel(cubepos_cen);
    for(int i = 0; i < 4; i++) cubechan_bei[i] = GetCubeChannel(cubepos_bei[i]);
 
    if(ADC_temp[int(cubechan_cen.Y())-1]<fADCCut || ADC_temp[int(cubechan_cen.X())-1]<fADCCut) continue;

    // Check which one used for noise estimation
    for(int i = 0; i < 4; i++){
      // Not include cube outside the matrix
      if(cubechan_bei[i].X()==-1){
        cubely_bei[i] = 0;
        continue;
      }

      // The left one would be used
      if(cubechan_bei[i].X()==cubechan_cen.X()) noise_temp = ADC_temp[int(cubechan_bei[i].Y())-1];
      else noise_temp = ADC_temp[int(cubechan_bei[i].X())-1];

      noise_esti->Fill(noise_temp);
      noise_ly->Fill(noise_temp);
    }*/

    // Method 2: estimate noise level per channel
    // Loop over each channel, if the channel satisfies requirement, then fill the noise histogram
    //if(fChanMapChoice==1){
      for(int i = 0; i < fChanNum; i++){
   
        // Channel ADC smaller than cut
        if(ADC_temp[i]>ADC_cut) continue;

        // Get the channel number of nearby channels
        channel_near = GetNearbyChannel(chan_order[i]);

        // All nearby channel ADCs should be smaller than cut value
        bool nearby_cut = true;
        if(channel_near.X()!=-1){
          if(data->ADC(channel_near.X()-1)>ADC_cut) nearby_cut = false;
        }
        if(channel_near.Y()!=-1){
          if(data->ADC(channel_near.Y()-1)>ADC_cut) nearby_cut = false;
        }
        if(channel_near.Z()!=-1){
          if(data->ADC(channel_near.Z()-1)>ADC_cut) nearby_cut = false;
        }
        if(channel_near.T()!=-1){
          if(data->ADC(channel_near.T()-1)>ADC_cut) nearby_cut = false;
        }

        if(nearby_cut==false) continue;

        // If all requirements are satisfied, fill the histogram
        noise_esti[i]->Fill(ADC_temp[i]);
        noise_ly->Fill(ADC_temp[i]);

      } // End of channel loop
    //}
    
    // Method 3: estimate noise level per channel, but in a safer way
    // This applies to the new board data, because there is trigger to select the events
    /*else if(fChanMapChoice==2){    
      // First match cubes at each layer, do not apply any cut
      cube_array = Get3DMatchedCubes(ADC_temp,ADC_cut);
      if(cube_array.size()!=3) continue;

      // Only select the vertical tracks
      //std::tie(track_info,chan_path,trk_graph) = Get3DTrackFit(cube_array);
      //if(TMath::Abs(track_info[1])>1e-4) continue;   
    
      if(CheckTrackVertical(cube_array)!=true) continue;
    
      // Loop over each layer
      for(int layer = 0; layer < 3; layer++){
        // Compute the distance of the cube to the center
        double dis = pow(cube_array[layer].X()-1,2) + pow(cube_array[layer].Y()-1,2);
        if(dis<1) continue;
      
        cubepos_cen.SetXYZ(cube_array[layer].X(),cube_array[layer].Y(),cube_array[layer].Z());
        cubepos_bei[0].SetXYZ(cube_array[layer].X()-2,cube_array[layer].Y(),cube_array[layer].Z());
        cubepos_bei[1].SetXYZ(cube_array[layer].X()+2,cube_array[layer].Y(),cube_array[layer].Z());
        cubepos_bei[2].SetXYZ(cube_array[layer].X(),cube_array[layer].Y()+2,cube_array[layer].Z());
        cubepos_bei[3].SetXYZ(cube_array[layer].X(),cube_array[layer].Y()-2,cube_array[layer].Z()); 
        for(int i = 0; i < 4; i++){
          cubechan_bei[i] = GetCubeChannel(cubepos_bei[i]);
          if(cubechan_bei[i].X()==-1) continue;
          // Fill the noise histogram, make sure the filled channel is not the same as track channel
          if(cubepos_bei[i].X()==cubepos_cen.X()){
            index_temp = ReturnIndex(int(cubechan_bei[i].Y()));
            noise_esti[index_temp]->Fill(data->ADC(int(cubechan_bei[i].Y())-1));
            noise_ly->Fill(data->ADC(int(cubechan_bei[i].Y())-1));
          }
          else if(cubepos_bei[i].Y()==cubepos_cen.Y()){
            index_temp = ReturnIndex(int(cubechan_bei[i].X()));
            noise_esti[index_temp]->Fill(data->ADC(int(cubechan_bei[i].X())-1));
            noise_ly->Fill(data->ADC(int(cubechan_bei[i].X())-1));
          }
        }
      }
    }*/

  }

  // Fit the noise value distribution with a Gaussian function
  /*double mean_temp = noise_esti->GetMean();
  double rms_temp = noise_esti->GetRMS();
  noise_esti->Fit("gaus","","",mean_temp-rms_temp,mean_temp+rms_temp);
  TF1 *fit_func = noise_esti->GetFunction("gaus");
  double noise_mean = fit_func->GetParameter(1);
  double noise_sigma = fit_func->GetParameter(2);*/

  // Noise distribution fit range (only consider the first peak)
  // Currently only used for SFGD new board data
  double noise_fitlow[fChanNum] = {100,90,83,73,77,54,87,30,83,73,92,83,54,80,79,72,89,79};
  double noise_fitupp[fChanNum] = {154,131,135,129,122,125,143,91,156,136,140,150,125,131,153,140,140,146};

  double mean_temp, rms_temp;
  double noise_mean[fChanNum];
  double noise_rms[fChanNum];
  TF1 *fit_func[fChanNum];
  for(int i = 0; i < fChanNum; i++){

    // First check whether the histogram is empty
    if(!((noise_esti[i]->Integral())>0)){
      noise_mean[i] = 0;
      noise_rms[i] = 0;
      continue;
    }

    // Get overall mean and RMS
    mean_temp = noise_esti[i]->GetMean();
    rms_temp = noise_esti[i]->GetRMS();

    //if(file_option==3){ // Gauss fit around the peak
    //  noise_esti[i]->Fit("gaus","","",noise_fitlow[i],noise_fitupp[i]);
    //  fit_func[i] = noise_esti[i]->GetFunction("gaus");
    //  noise_mean[i] = fit_func[i]->GetParameter(1);
    //  noise_rms[i] = fit_func[i]->GetParameter(2);
    //}
    //else{ // Take the overall mean and RMS
      noise_mean[i] = mean_temp;
      noise_rms[i] = rms_temp; 
    //}
   
  }
  std::cout << "Noise level estimation done" << std::endl;

  // Estimate the crosstalk fraction
  std::cout << "Start to estimate crosstalk rate" << std::endl;
  // Second loop over all events 
  for(int n = 0; n < n_event; n++){

    data->GetMppc(n,swap,newswap);

    ADC_temp.clear();
    for(int i = 0; i < fChanNum; i++) ADC_temp.push_back(data->ADC(chan_order[i]-1));

    // Pass cosmic cut
    pass_tag = EventCosmicCut(ADC_temp,ADC_cut);
    if(pass_tag==false) continue;
 
    cube_array = Get3DMatchedCubes(ADC_temp,ADC_cut);
    if(cube_array.size()!=3) continue;

    // Only select the vertical tracks
    //std::tie(track_info,chan_path,trk_graph) = Get3DTrackFit(cube_array);
    //if(TMath::Abs(track_info[1])>1e-4) continue;
  
    if(CheckTrackVertical(cube_array)!=true) continue;

    // Only select the edge tracks
    double dis = pow(cube_array[1].X()-1,2) + pow(cube_array[1].Y()-1,2);
    if(dis!=1.) continue;

    // Estimate the crosstalk rate (loop over all 3 layers)
    for(int layer = 0; layer < 3; layer++){

    cubepos_cen.SetXYZ(cube_array[layer].X(),cube_array[layer].Y(),cube_array[layer].Z());
    cubepos_bei[0].SetXYZ(cube_array[layer].X()-1,cube_array[layer].Y(),cube_array[layer].Z());
    cubepos_bei[1].SetXYZ(cube_array[layer].X()+1,cube_array[layer].Y(),cube_array[layer].Z());
    cubepos_bei[2].SetXYZ(cube_array[layer].X(),cube_array[layer].Y()+1,cube_array[layer].Z());
    cubepos_bei[3].SetXYZ(cube_array[layer].X(),cube_array[layer].Y()-1,cube_array[layer].Z());

    cubechan_cen = GetCubeChannel(cubepos_cen);
    for(int i = 0; i < 4; i++) cubechan_bei[i] = GetCubeChannel(cubepos_bei[i]);

    // Do a ADC cut, should not apply on the middle layer
    if(data->ADC(int(cubechan_cen.Y())-1)<ADC_cut && layer!=1) continue;
    if(data->ADC(int(cubechan_cen.X())-1)<ADC_cut && layer!=1) continue;

    // Check how many nearby cubes in each side (X direction or Y direction)
    nx = 0; ny = 0;
    for(int i = 0; i < 4; i++){
      // Not include cube outside the matrix
      if(cubechan_bei[i].X()==-1) continue;

      // Check X direction
      if(cubepos_bei[i].X()==cubepos_cen.X()) nx += 1;

      // Check Y direction
      if(cubepos_bei[i].Y()==cubepos_cen.Y()) ny += 1;
    }

    // Method 1: consider cube LY from two fibers
    // First estimate the light yield in each nearby cube
    /*for(int i = 0; i < 4; i++){
      // Not include cube outside the matrix
      if(cubechan_bei[i].X()==-1){
        cubely_bei[i] = 0;
        continue;
      } 

      if(cubechan_bei[i].X()==cubechan_cen.X()) cubely_bei[i] = ADC_temp[int(cubechan_bei[i].Y())-1] * 2;
      else cubely_bei[i] = ADC_temp[int(cubechan_bei[i].X())-1] * 2;
    }

    // Then estimate the light yield in the center cube
    cubely_cen = (ADC_temp[int(cubechan_cen.X())-1] - (cubely_bei[2] + cubely_bei[3]) / 2)
                 + (ADC_temp[int(cubechan_cen.Y())-1] - (cubely_bei[0] + cubely_bei[1]) / 2);

    // Fill the histograms
    for(int i = 0; i < 4; i++){
      // Not include cube outside the matrix
      if(cubechan_bei[i].X()==-1) continue;

      xtalk_rate->Fill(cubely_bei[i]/cubely_cen*100);
      beicube_ly->Fill(cubely_bei[i]);
    }
    trkcube_ly->Fill(cubely_cen);*/

    // Method 2: consider cube LY from only one fiber
    for(int i = 0; i < 4; i++){
      // Not include cube outside the matrix
      if(cubechan_bei[i].X()==-1){
        cubely_bei[i] = 0;
        cubenoise_bei[i] = 0;
        cubegain_bei[i] = 0;
        continue;
      } 

      if(cubechan_bei[i].X()==cubechan_cen.X()){
        cubely_bei[i] = data->ADC(int(cubechan_bei[i].Y())-1);
        index_temp = ReturnIndex(int(cubechan_bei[i].Y()));
        cubenoise_bei[i] = noise_mean[index_temp]; // Estimate noise level per channel
        cubegain_bei[i] = ADCgain[index_temp]; // Use gain per channel     
        xtalk_adc[index_temp]->Fill(data->ADC(int(cubechan_bei[i].Y())-1));   
      }
      else if(cubechan_bei[i].Y()==cubechan_cen.Y()){      
        cubely_bei[i] = data->ADC(int(cubechan_bei[i].X())-1);
        index_temp = ReturnIndex(int(cubechan_bei[i].X()));
        cubenoise_bei[i] = noise_mean[index_temp]; // Estimate noise level per channel
        cubegain_bei[i] = ADCgain[index_temp]; // Use gain per channel
        xtalk_adc[index_temp]->Fill(data->ADC(int(cubechan_bei[i].X())-1));
      }
    }

    // Compute crosstalk
    for(int i = 0; i < 4; i++){
      // Not include cube outside the matrix
      if(cubechan_bei[i].X()==-1) continue;

      if(cubechan_bei[i].X()==cubechan_cen.X()){
        if(data->ADC(int(cubechan_cen.Y())-1)>fADCUppCut) continue;
        if(cubenoise_bei[i]<=0.) continue; 

        //cubely_cen = ADC_temp[int(cubechan_cen.Y())-1] - (cubely_bei[0] + cubely_bei[1]);
        //xtalk_frac = (cubely_bei[i] - noise_mean) / (ADC_temp[int(cubechan_cen.Y())-1] - 2 * cubely_bei[i] + noise_mean) * 100;
        //xtalk_frac = (cubely_bei[i] - noise_mean) / (ADC_temp[int(cubechan_cen.Y())-1] - noise_mean) * 100;

        index_temp = ReturnIndex(int(cubechan_cen.Y()));
        xtalk_frac = (cubely_bei[i] - cubenoise_bei[i]) / (data->ADC(int(cubechan_cen.Y())-1) - noise_mean[index_temp]) * 100;
        noise_frac = cubenoise_bei[i] / (data->ADC(int(cubechan_cen.Y())-1) - noise_mean[index_temp]) * 100;

        // Crosstalk level per channel
        if(xtalk_frac>=0) xtalk_esti[index_temp]->Fill(xtalk_frac);
        else xtalk_esti[index_temp]->Fill(0);

        // Only use middle layer to estimate crosstalk
        if(layer==1){
          if(xtalk_frac>=0) xtalk_rate->Fill(xtalk_frac);
          else xtalk_rate->Fill(0);
          //xtalk_rate->Fill(xtalk_frac);
          comp_noisefrac->Fill(noise_frac);
          xtalk_frac = cubely_bei[i] / (data->ADC(int(cubechan_cen.Y())-1) - noise_mean[index_temp]) * 100;
          comp_xtalkfrac->Fill(xtalk_frac);
        }

        beicube_ly->Fill(cubely_bei[i]);
        trkcube_ly->Fill(data->ADC(int(cubechan_cen.Y())-1));
        track_adc[index_temp]->Fill(data->ADC(int(cubechan_cen.Y())-1));

        // Crosstalk level with p.e
        xtalk_frac = ((cubely_bei[i] - cubenoise_bei[i]) / cubegain_bei[i]) / ((data->ADC(int(cubechan_cen.Y())-1) - noise_mean[index_temp]) / ADCgain[index_temp]) * 100;
        noise_frac = (cubenoise_bei[i] / cubegain_bei[i]) / ((data->ADC(int(cubechan_cen.Y())-1) - noise_mean[index_temp]) / ADCgain[index_temp]) * 100;
        if(layer==1){
          double trk_pe = (data->ADC(int(cubechan_cen.Y())-1) - noise_mean[index_temp]) / ADCgain[index_temp];
          double xtalk_pe = (cubely_bei[i] - cubenoise_bei[i]) / cubegain_bei[i];
          if(xtalk_pe<=0) xtalk_pe = 0;
          test_2d->Fill(trk_pe,xtalk_pe);
          if(xtalk_frac>=0) xtalk_rate_pe->Fill(xtalk_frac);
          else xtalk_rate_pe->Fill(0);
          //xtalk_rate_pe->Fill(xtalk_frac);
          comp_noisefrac_pe->Fill(noise_frac);
          xtalk_frac = (cubely_bei[i] / cubegain_bei[i]) / ((data->ADC(int(cubechan_cen.Y())-1) - noise_mean[index_temp]) / ADCgain[index_temp]) * 100;
          comp_xtalkfrac_pe->Fill(xtalk_frac);
        }
      }
      else if(cubechan_bei[i].Y()==cubechan_cen.Y()){
        if(data->ADC(int(cubechan_cen.X())-1)>fADCUppCut) continue; 
        if(cubenoise_bei[i]<=0.) continue; 

        //cubely_cen = ADC_temp[int(cubechan_cen.X())-1] - (cubely_bei[2] + cubely_bei[3]);
        //xtalk_frac = (cubely_bei[i] - noise_mean) / (ADC_temp[int(cubechan_cen.X())-1] - 2 * cubely_bei[i] + noise_mean) * 100;
        //xtalk_frac = (cubely_bei[i] - noise_mean) / (ADC_temp[int(cubechan_cen.X())-1] - noise_mean) * 100;

        index_temp = ReturnIndex(int(cubechan_cen.X()));
        xtalk_frac = (cubely_bei[i] - cubenoise_bei[i]) / (data->ADC(int(cubechan_cen.X())-1) - noise_mean[index_temp]) * 100;
        noise_frac = cubenoise_bei[i] / (data->ADC(int(cubechan_cen.X())-1) - noise_mean[index_temp]) * 100;

        // Crosstalk level per channel
        if(xtalk_frac>=0) xtalk_esti[index_temp]->Fill(xtalk_frac);
        else xtalk_esti[index_temp]->Fill(0);

        // Only use middle layer to estimate crosstalk
        if(layer==1){ 
          if(xtalk_frac>=0) xtalk_rate->Fill(xtalk_frac);
          else xtalk_rate->Fill(0);
          //xtalk_rate->Fill(xtalk_frac);
          comp_noisefrac->Fill(noise_frac);
          xtalk_frac = cubely_bei[i] / (data->ADC(int(cubechan_cen.X())-1) - noise_mean[index_temp]) * 100;
          comp_xtalkfrac->Fill(xtalk_frac);
        }

        beicube_ly->Fill(cubely_bei[i]);
        trkcube_ly->Fill(data->ADC(int(cubechan_cen.X())-1));
        track_adc[index_temp]->Fill(data->ADC(int(cubechan_cen.X())-1));

        // Crosstalk level with p.e.
        xtalk_frac = ((cubely_bei[i] - cubenoise_bei[i]) / cubegain_bei[i]) / ((data->ADC(int(cubechan_cen.X())-1) - noise_mean[index_temp]) / ADCgain[index_temp]) * 100;
        noise_frac = (cubenoise_bei[i] / cubegain_bei[i]) / ((data->ADC(int(cubechan_cen.X())-1) - noise_mean[index_temp]) / ADCgain[index_temp]) * 100;
        if(layer==1){
          double trk_pe = (data->ADC(int(cubechan_cen.X())-1) - noise_mean[index_temp]) / ADCgain[index_temp];
          double xtalk_pe = (cubely_bei[i] - cubenoise_bei[i]) / cubegain_bei[i];
          if(xtalk_pe<=0) xtalk_pe = 0;
          test_2d->Fill(trk_pe,xtalk_pe);
          if(xtalk_frac>=0) xtalk_rate_pe->Fill(xtalk_frac);
          else xtalk_rate_pe->Fill(0);
          //xtalk_rate_pe->Fill(xtalk_frac);
          comp_noisefrac_pe->Fill(noise_frac);
          xtalk_frac = (cubely_bei[i] / cubegain_bei[i]) / ((data->ADC(int(cubechan_cen.X())-1) - noise_mean[index_temp]) / ADCgain[index_temp]) * 100;
          comp_xtalkfrac_pe->Fill(xtalk_frac);
        }
      }
    }

    } // End of three layer loop

  }

  std::cout << "Crosstalk rate estimation done" << std::endl;

  // Fit the crosstalk distribution with Landau function
  double range_low = xtalk_rate->GetMean() - 2 * xtalk_rate->GetRMS();
  double range_upp = xtalk_rate->GetMean() + 4 * xtalk_rate->GetRMS();
  xtalk_rate->Fit("landau","","",range_low,range_upp);
  TF1 *fitfunc_xtalk = xtalk_rate->GetFunction("landau");
  double mean = fitfunc_xtalk->GetParameter(1);
  double rms = fitfunc_xtalk->GetParameter(2);

  range_low = xtalk_rate_pe->GetMean() - 2 * xtalk_rate_pe->GetRMS();
  range_upp = xtalk_rate_pe->GetMean() + 4 * xtalk_rate_pe->GetRMS();
  xtalk_rate_pe->Fit("landau","","",range_low,range_upp);
  fitfunc_xtalk = xtalk_rate_pe->GetFunction("landau");
  double mean_ex = fitfunc_xtalk->GetParameter(1);
  double rms_ex = fitfunc_xtalk->GetParameter(2);

  std::cout << "Crosstalk rate MPV (ADC): " << mean << endl;
  std::cout << "Crosstalk rate MPV (p.e.): " << mean_ex << endl;

  TText *pl_name = new TText();
  pl_name->SetTextSize(0.04);
  TText *pl_mean = new TText();
  pl_mean->SetTextSize(0.03);
  TText *pl_rms = new TText();
  pl_rms->SetTextSize(0.03);

  // Scale the histograms
  double scale_trkcube = trkcube_ly->Integral();
  double scale_beicube = beicube_ly->Integral();
  double scale_noise = noise_ly->Integral();
  if(scale_trkcube>0) trkcube_ly->Scale(1./scale_trkcube);
  if(scale_beicube>0) beicube_ly->Scale(1./scale_beicube);
  if(scale_noise>0) noise_ly->Scale(1./scale_noise);
  trkcube_ly->GetYaxis()->SetRangeUser(1e-3,1);
  beicube_ly->GetYaxis()->SetRangeUser(1e-3,1);
  noise_ly->GetYaxis()->SetRangeUser(1e-3,1);

  // Print the noise level
  std::cout << "Noise level: {";
  for(int i = 0; i < fChanNum; i++) std::cout << noise_mean[i] << ", ";
  std::cout << "}" << endl;

  // Draw the plots
  gStyle->SetOptStat(0);

  TCanvas *c0 = new TCanvas("frac_comp","frac_comp",700,600);
  c0->SetLeftMargin(0.15);
  c0->cd();
  comp_noisefrac->SetLineWidth(2);
  comp_noisefrac->SetLineColor(kBlue);
  comp_xtalkfrac->SetLineWidth(2);
  comp_xtalkfrac->SetLineColor(kRed);
  comp_noisefrac->Draw("hist");
  comp_xtalkfrac->Draw("hist same");
  
  TLegend *lg0 = new TLegend(0.6,0.68,0.89,0.87);
  lg0->SetLineWidth(0);
  lg0->SetFillStyle(0);
  lg0->AddEntry(comp_noisefrac,"Noise fraction","l");
  lg0->AddEntry(comp_xtalkfrac,"Crosstalk fraction","l");
  lg0->Draw("same");
  
  gPad->SetGridx();
  gPad->SetGridy();
  c0->Update(); 
  
  TCanvas *c0ex = new TCanvas("frac_comp_pe","frac_comp_pe",700,600);
  c0ex->SetLeftMargin(0.15);
  c0ex->cd();
  comp_noisefrac_pe->SetLineWidth(2);
  comp_noisefrac_pe->SetLineColor(kBlue);
  comp_xtalkfrac_pe->SetLineWidth(2);
  comp_xtalkfrac_pe->SetLineColor(kRed);
  comp_noisefrac_pe->Draw("hist");
  comp_xtalkfrac_pe->Draw("hist same");
  
  TLegend *lg0ex = new TLegend(0.6,0.68,0.89,0.87);
  lg0ex->SetLineWidth(0);
  lg0ex->SetFillStyle(0);
  lg0ex->AddEntry(comp_noisefrac_pe,"Noise fraction","l");
  lg0ex->AddEntry(comp_xtalkfrac_pe,"Crosstalk fraction","l");
  lg0ex->Draw("same");
  
  gPad->SetGridx();
  gPad->SetGridy();
  c0ex->Update(); 

  TCanvas *c1 = new TCanvas("xtalk_rate","xtalk_rate",700,600);
  c1->SetLeftMargin(0.15);
  c1->cd();
  xtalk_rate->SetLineWidth(2);
  xtalk_rate->SetLineColor(kBlue);
  //xtalk_rate->SetMarkerColor(kBlue);
  //xtalk_rate->SetMarkerStyle(20);
  //xtalk_rate->Draw("P E0");
  xtalk_rate->Draw("hist");

  name.Form("Crosstalk fraction:");
  pl_name->DrawTextNDC(0.5,0.68,name);
  name.Form("Mean = %f",xtalk_rate->GetMean());
  pl_mean->DrawTextNDC(0.5,0.61,name);
  gPad->SetGridx();
  gPad->SetGridy();
  if(file_option==1 || file_option==4) gPad->SetLogy();
  c1->Update();

  TCanvas *c1ex = new TCanvas("xtalk_rate_pe","xtalk_rate_pe",700,600);
  c1ex->SetLeftMargin(0.15);
  c1ex->cd();
  xtalk_rate_pe->SetLineWidth(2);
  xtalk_rate_pe->SetLineColor(kRed);
  //xtalk_rate_pe->SetMarkerColor(kRed);
  //xtalk_rate_pe->SetMarkerStyle(20);
  //xtalk_rate_pe->Draw("E P0");
  xtalk_rate_pe->Draw("hist");

  name.Form("Crosstalk fraction:");
  pl_name->DrawTextNDC(0.5,0.68,name);
  name.Form("Mean = %f",xtalk_rate_pe->GetMean());
  pl_mean->DrawTextNDC(0.5,0.61,name);
  gPad->SetGridx();
  gPad->SetGridy();
  if(file_option==1 || file_option==4) gPad->SetLogy();
  c1ex->Update();

  TCanvas *c2 = new TCanvas("cubely_comp","cubely_comp",700,600);
  c2->SetLeftMargin(0.15);
  c2->cd();

  trkcube_ly->SetLineWidth(2);
  trkcube_ly->SetLineColor(kRed);
  beicube_ly->SetLineWidth(2);
  beicube_ly->SetLineColor(kBlue);
  noise_ly->SetLineWidth(2);
  noise_ly->SetLineColor(kOrange-7);
  
  trkcube_ly->Draw("hist");
  beicube_ly->Draw("hist same");
  noise_ly->Draw("hist same");
  TLegend *lg2 = new TLegend(0.35,0.6,0.7,0.87);
  lg2->SetLineWidth(0);
  lg2->SetFillStyle(0);
  lg2->AddEntry(trkcube_ly,"Track channel","l");
  lg2->AddEntry(beicube_ly,"Crosstalk channel","l");
  lg2->AddEntry(noise_ly,"Noise level","l");
  lg2->Draw("same");
  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();
  c2->Update();

  /*TCanvas *c3 = new TCanvas("noise_esti","noise_esti",700,600);
  c3->SetLeftMargin(0.15);
  c3->cd();
  noise_esti->SetLineWidth(2);
  noise_esti->SetLineColor(kBlue);
  noise_esti->Draw("hist");

  name.Form("Gauss fit:");
  pl_name->DrawTextNDC(0.55,0.68,name);
  name.Form("Mean = %f ADC",noise_mean);
  pl_mean->DrawTextNDC(0.55,0.61,name);
  gPad->SetGridx();
  gPad->SetGridy();
  c3->Update();*/

  // Drawing order array
  // Noise level
  double st_left = 0.55;
  double st_top = 0.7;
  pl_mean->SetTextSize(0.04);
  pl_rms->SetTextSize(0.04);
 
  // MPPC Face 1 + 3
  TCanvas *c4 = new TCanvas("noise_xz","noise_xz",1200,1200);
  c4->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c4->cd(i+1);
    noise_esti[i]->Draw("hist");
    //if(file_option==3 && noise_mean[i]!=0) fit_func[i]->Draw("same");
    name.Form("Overall mean = %f ADC",noise_mean[i]);
    pl_mean->DrawTextNDC(st_left,st_top,name);
    name.Form("Overall RMS = %f ADC",noise_rms[i]);
    pl_rms->DrawTextNDC(st_left,st_top-0.07,name);
    //gPad->SetLogy();
  }
  c4->Update();

  // MPPC Face 2 + 4
  TCanvas *c5 = new TCanvas("noise_yz","noise_yz",1200,1200);
  c5->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c5->cd(i+1);
    noise_esti[i+9]->Draw("hist");
    //if(file_option==3 && noise_mean[i]!=0) fit_func[i+9]->Draw("same");
    name.Form("Overall mean = %f ADC",noise_mean[i+9]);
    pl_mean->DrawTextNDC(st_left,st_top,name);
    name.Form("Overall RMS = %f ADC",noise_rms[i+9]);
    pl_rms->DrawTextNDC(st_left,st_top-0.07,name);
    //gPad->SetLogy();
  }
  c5->Update();

  // Crosstalk level
  // Face 1 + 3
  TCanvas *c6 = new TCanvas("xtalk_xz","xtalk_xz",1200,1200);
  c6->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c6->cd(i+1);
    xtalk_esti[i]->Draw("hist");
  }
  c6->Update();
  
  // Face 2 + 4
  TCanvas *c7 = new TCanvas("xtalk_yz","xtalk_yz",1200,1200);
  c7->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c7->cd(i+1);
    xtalk_esti[i+9]->Draw("hist");
  }
  c7->Update();

  // Crosstalk channel ADC
  // Face 1 + 3
  TCanvas *c8 = new TCanvas("xtalkADC_xz","xtalkADC_xz",1200,1200);
  c8->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c8->cd(i+1);
    xtalk_adc[i]->Draw("hist");
  }
  c8->Update();

  // Face 2 + 4
  TCanvas *c9 = new TCanvas("xtalkADC_yz","xtalkADC_yz",1200,1200);
  c9->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c9->cd(i+1);
    xtalk_adc[i+9]->Draw("hist");
  }
  c9->Update();

  // Track channel ADC
  // Face 1 + 3
  TCanvas *c10 = new TCanvas("trackADC_xz","trackADC_xz",1200,1200);
  c10->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c10->cd(i+1);
    track_adc[i]->Draw("hist");
  }
  c10->Update();

  // Face 2 + 4
  TCanvas *c11 = new TCanvas("trackADC_yz","trackADC_yz",1200,1200);
  c11->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c11->cd(i+1);
    track_adc[i+9]->Draw("hist");
  }
  c11->Update();

  TCanvas *c12 = new TCanvas("test_2d","test_2d",800,600);
  c12->SetLeftMargin(0.15);
  c12->SetRightMargin(0.15);
  c12->cd();
  test_2d->SetContour(99);
  test_2d->Draw("colz");
  c12->Update();

  TString prefix = "../../../plots/scintillator_cube/";
  TString type = fPathName[file_option-1] + "/";
  TString suffix;

  suffix = prefix + type + "frac_comp.png";
  c0->SaveAs(suffix);

  suffix = prefix + type + "frac_comp_pe.png";
  c0ex->SaveAs(suffix);

  suffix = prefix + type + "xtalk_rate.png";
  c1->SaveAs(suffix);

  suffix = prefix + type + "xtalk_rate.pdf";
  c1->SaveAs(suffix);

  suffix = prefix + type + "xtalk_rate_pe.png";
  c1ex->SaveAs(suffix);

  suffix = prefix + type + "xtalk_rate_pe.pdf";
  c1ex->SaveAs(suffix);

  suffix = prefix + type + "cubely_comp.png";
  c2->SaveAs(suffix);

  //suffix = prefix + type + "noise_esti.png";
  //c3->SaveAs(suffix);

  suffix = prefix + type + "noise_xz.png";
  c4->SaveAs(suffix);

  suffix = prefix + type + "noise_yz.png";
  c5->SaveAs(suffix);

  suffix = prefix + type + "xtalk_xz.png";
  c6->SaveAs(suffix);

  suffix = prefix + type + "xtalk_yz.png";
  c7->SaveAs(suffix);

  suffix = prefix + type + "xtalkADC_xz.png";
  c8->SaveAs(suffix);

  suffix = prefix + type + "xtalkADC_yz.png";
  c9->SaveAs(suffix);

  suffix = prefix + type + "trackADC_xz.png";
  c10->SaveAs(suffix);

  suffix = prefix + type + "trackADC_yz.png";
  c11->SaveAs(suffix);

  suffix = prefix + type + "test_2d.png";
  c12->SaveAs(suffix);

  // Save the plots into output file
  TString fout_name = "../../results/CrosstalkAnalysis_" + fOutFileName[file_option-1] + ".root";

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();
 
  c0->Write();
  c0ex->Write();
  c1->Write();
  c1ex->Write();
  c2->Write();
  //c3->Write();
  c4->Write();
  c5->Write();
  c6->Write();
  c7->Write();
  c8->Write();
  c9->Write();
  c10->Write();
  c11->Write();
  c12->Write();

  xtalk_rate->Write();
  xtalk_rate_pe->Write();
  trkcube_ly->Write();
  beicube_ly->Write();
  //noise_esti->Write();

  for(int i = 0; i < fChanNum; i++){
    noise_esti[i]->Write();
    xtalk_esti[i]->Write();
    xtalk_adc[i]->Write();
    track_adc[i]->Write();
  }

  test_2d->Write();

  fout->Close();

}

// General event 3D analysis
void Event3DAnalysis(int file_option = 1, double ADC_cut = fADCCut){

  std::string fin_name;
  bool swap, newswap;

  // Set input file name
  fin_name = "../../inputs/" + fFileName[file_option-1] + ".root";
  swap = fSwap[file_option-1];
  newswap = fNewSwap[file_option-1];

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();
  
  std::cout << "Total number of events: " << n_event << std::endl;
  
  // Channel mapping choice
  if(file_option==1 || file_option==2) fChanMapChoice = 1;
  else if(file_option==3 || file_option==4) fChanMapChoice = 2;

  int chan_order[fChanNum];
  for(int i = 0; i < fChanNum; i++) chan_order[i] = fChanOrder[fChanMapChoice-1][i];

  // Create some variables
  int n_pass = 0;
  bool pass_tag;
  std::vector<double> ADC_temp;
  std::vector<TLorentzVector> cube_array;
  double tot_ly;
  double trk_len;
  double trk_ang_pol, trk_ang_azi;
  TGraph2D *trk_graph;
  std::vector<double> track_info;
  std::vector<double> chan_path;

  // Create some histograms
  double totly_lim = 30000;

  // Total light yield vs. track direction (polar angle)
  TH2D *totly_ang_polar = new TH2D("totly_ang_polar","totly_ang_polar",10,0,45,60,0,totly_lim);
  totly_ang_polar->GetXaxis()->SetTitle("Track polar angle #theta / degree");
  totly_ang_polar->GetYaxis()->SetTitle("Track light yield / ADC");
  totly_ang_polar->GetXaxis()->SetLabelSize(0.04);
  totly_ang_polar->GetXaxis()->SetTitleSize(0.04);
  totly_ang_polar->GetYaxis()->SetLabelSize(0.04);
  totly_ang_polar->GetYaxis()->SetTitleSize(0.04);
  totly_ang_polar->GetYaxis()->SetTitleOffset(1.4);
  totly_ang_polar->GetZaxis()->SetLabelSize(0.04);
  totly_ang_polar->SetTitle("");

  // (Azimuth angle)
  TH2D *totly_ang_azi = new TH2D("totly_ang_azi","totly_ang_azi",20,-180,180,60,0,totly_lim);
  totly_ang_azi->GetXaxis()->SetTitle("Track azimuth angle #phi / degree");
  totly_ang_azi->GetYaxis()->SetTitle("Track light yield / ADC");
  totly_ang_azi->GetXaxis()->SetLabelSize(0.04);
  totly_ang_azi->GetXaxis()->SetTitleSize(0.04);
  totly_ang_azi->GetYaxis()->SetLabelSize(0.04);
  totly_ang_azi->GetYaxis()->SetTitleSize(0.04);
  totly_ang_azi->GetYaxis()->SetTitleOffset(1.4);
  totly_ang_azi->GetZaxis()->SetLabelSize(0.04);
  totly_ang_azi->SetTitle("");

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

  // Distribution of angles (polar + azimuth)
  TH1D *ang_polar = new TH1D("ang_polar","ang_polar",10,0,45);
  ang_polar->GetXaxis()->SetTitle("Track polar angle #theta / degree");
  ang_polar->GetYaxis()->SetTitle("Number of events / bin");
  ang_polar->GetXaxis()->SetLabelSize(0.04);
  ang_polar->GetXaxis()->SetTitleSize(0.04);
  ang_polar->GetYaxis()->SetLabelSize(0.04);
  ang_polar->GetYaxis()->SetTitleSize(0.04);
  ang_polar->GetYaxis()->SetTitleOffset(1.4);
  ang_polar->SetTitle("");
 
  TH1D *ang_azi = new TH1D("ang_azi","ang_azi",20,-180,180);
  ang_azi->GetXaxis()->SetTitle("Track azimuth angle #phi / degree");
  ang_azi->GetYaxis()->SetTitle("Number of events / bin");
  ang_azi->GetXaxis()->SetLabelSize(0.04);
  ang_azi->GetXaxis()->SetTitleSize(0.04);
  ang_azi->GetYaxis()->SetLabelSize(0.04);
  ang_azi->GetYaxis()->SetTitleSize(0.04);
  ang_azi->GetYaxis()->SetTitleOffset(1.4);
  ang_azi->SetTitle("");

  // Distribution of distance sum (check track fit quality)
  TH1D *trkfit_quality = new TH1D("trkfit_quality","trkfit_quality",80,0,30);
  trkfit_quality->SetTitle("");
  trkfit_quality->GetXaxis()->SetTitle("Distance sum / mm");
  trkfit_quality->GetYaxis()->SetTitle("Number of events / bin");
  trkfit_quality->GetXaxis()->SetLabelSize(0.04);
  trkfit_quality->GetXaxis()->SetTitleSize(0.04);
  trkfit_quality->GetYaxis()->SetLabelSize(0.04);
  trkfit_quality->GetYaxis()->SetTitleSize(0.04);
  trkfit_quality->GetYaxis()->SetTitleOffset(1.4);
  trkfit_quality->SetTitle("");

  // Distribution of path length seen by each channel
  TH1D *path_length[fChanNum];
  TString title;
  for(int i = 0; i < fChanNum; i++){
    title.Form("channel%d_path",i+1);
    path_length[i] = new TH1D(title,title,60,0,20);
    path_length[i]->GetXaxis()->SetTitle("Estimated path length / mm");
    path_length[i]->GetYaxis()->SetTitle("Number of events / bin");
    path_length[i]->SetTitle(title);
    path_length[i]->SetLineWidth(2);
    path_length[i]->SetLineColor(kBlue);
  }

  // Loop over all events
  for(int n = 0; n < n_event; n++){

    data->GetMppc(n,swap,newswap);

    ADC_temp.clear();
    for(int i = 0; i < fChanNum; i++) ADC_temp.push_back(data->ADC(chan_order[i]-1));

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
    std::tie(track_info,chan_path,trk_graph) = Get3DTrackFit(cube_array);
    trk_len = track_info[0];
    trk_ang_pol = track_info[1];
    trk_ang_azi = track_info[2];

    // Fill the histograms
    totly_ang_polar->Fill(trk_ang_pol,tot_ly);
    totly_ang_azi->Fill(trk_ang_azi,tot_ly);
    locally->Fill(tot_ly/trk_len);
    ang_polar->Fill(trk_ang_pol);
    ang_azi->Fill(trk_ang_azi);

    trkfit_quality->Fill(track_info[3]);
    for(int i = 0; i < fChanNum; i++) path_length[i]->Fill(chan_path[i]);

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

  std::cout << "Local light yield MPV: " << mean << endl;

  TString name;
  TText *pl_name = new TText();
  pl_name->SetTextSize(0.04);
  TText *pl_mean = new TText();
  pl_mean->SetTextSize(0.03);
  TText *pl_rms = new TText();
  pl_rms->SetTextSize(0.03);
 
  // Draw the plots
  gStyle->SetOptStat(0);
  double left_begin = 0.55;
  double top_begin = 0.68;

  TCanvas *c1 = new TCanvas("totly_ang_polar","totly_ang_polar",800,600);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.15);
  c1->cd();
  totly_ang_polar->SetContour(99);
  totly_ang_polar->Draw("colz");
  c1->Update();

  TCanvas *c2 = new TCanvas("totly_ang_azi","totly_ang_azi",800,600);
  c2->SetLeftMargin(0.15);
  c2->SetRightMargin(0.15);
  c2->cd();
  totly_ang_azi->SetContour(99);
  totly_ang_azi->Draw("colz");
  c2->Update();

  TCanvas *c3 = new TCanvas("locally","locally",700,600);
  c3->SetLeftMargin(0.15);
  c3->cd();
  locally->SetLineWidth(2);
  locally->SetLineColor(kBlue);
  locally->Draw("hist");

  name.Form("Landau fit:");
  pl_name->DrawTextNDC(left_begin,top_begin,name);
  name.Form("MPV = %f ADC / mm",mean);
  pl_mean->DrawTextNDC(left_begin,top_begin-0.07,name);

  gPad->SetGridx();
  gPad->SetGridy();
  c3->Update();

  TCanvas *c4 = new TCanvas("ang_polar","ang_polar",700,600);
  c4->SetLeftMargin(0.15);
  c4->cd();
  ang_polar->SetLineWidth(2);
  ang_polar->SetLineColor(kBlue);
  ang_polar->Draw("hist");
  gPad->SetGridx();
  gPad->SetGridy();
  c4->Update();

  TCanvas *c5 = new TCanvas("ang_azi","ang_azi",700,600);
  c5->SetLeftMargin(0.15);
  c5->cd();
  ang_azi->SetLineWidth(2);
  ang_azi->SetLineColor(kBlue);
  ang_azi->Draw("hist");
  gPad->SetGridx();
  gPad->SetGridy();
  c5->Update();

  TCanvas *c6 = new TCanvas("trkfit_quality","trkfit_quality",700,600);
  c6->SetLeftMargin(0.15);
  c6->cd();
  trkfit_quality->SetLineWidth(2);
  trkfit_quality->SetLineColor(kBlue);
  trkfit_quality->Draw("hist");
  gPad->SetGridx();
  gPad->SetGridy();
  c6->Update();

  // MPPC channel face 1 + 3
  TCanvas *c7 = new TCanvas("pathlength_xz","pathlength_xz",1200,1200);
  c7->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c7->cd(i+1);
    path_length[i]->Draw("hist");
  }
  c7->Update();

  // MPPC channel face 2 + 4
  TCanvas *c8 = new TCanvas("pathlength_yz","pathlength_yz",1200,1200);
  c8->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c8->cd(i+1);
    path_length[i+9]->Draw("hist");
  }
  c8->Update();
  
  TString prefix = "../../../plots/scintillator_cube/";
  TString type = fPathName[file_option-1] + "/";
  TString suffix;

  suffix = prefix + type + "totly_ang_polar.png";
  c1->SaveAs(suffix);

  suffix = prefix + type + "totly_ang_azi.png";
  c2->SaveAs(suffix);

  suffix = prefix + type + "locally.png";
  c3->SaveAs(suffix);

  suffix = prefix + type + "ang_polar.png";
  c4->SaveAs(suffix);

  suffix = prefix + type + "ang_azi.png";
  c5->SaveAs(suffix);

  suffix = prefix + type + "trkfit_quality.png";
  c6->SaveAs(suffix);

  suffix = prefix + type + "pathlength_xz.png";
  c7->SaveAs(suffix);

  suffix = prefix + type + "pathlength_yz.png";
  c8->SaveAs(suffix);

  // Save the plots into output file
  TString fout_name = "../../results/Track3DAnalysis_" + fOutFileName[file_option-1] + ".root";

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

  totly_ang_polar->Write();
  totly_ang_azi->Write();
  locally->Write();
  ang_polar->Write();
  ang_azi->Write();

  trkfit_quality->Write();
  for(int i = 0; i < fChanNum; i++) path_length[i]->Write();

  fout->Close();

}

// Draw event 3D display (after some cuts)
// Notation: file_option = input data file want to be used
// Notation: seed = random number generator
// Notation: ADC_cut = a cut applied on ADC value
void DrawEvent3D(int file_option = 1, int seed = 0, double ADC_cut = fADCCut){

  std::string fin_name;
  bool swap, newswap;

  // Set input file name
  fin_name = "../../inputs/" + fFileName[file_option-1] + ".root";
  swap = fSwap[file_option-1];
  newswap = fNewSwap[file_option-1];

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();
  std::cout << "Total number of events: " << n_event << std::endl;
  int rand_event;

  // Channel mapping choice
  if(file_option==1 || file_option==2) fChanMapChoice = 1;
  else if(file_option==3 || file_option==4) fChanMapChoice = 2;

  int chan_order[fChanNum];
  for(int i = 0; i < fChanNum; i++) chan_order[i] = fChanOrder[fChanMapChoice-1][i];

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
  double trk_ang_pol, trk_ang_azi;
  std::vector<double> track_info;
  std::vector<double> chan_path;

  // Read the events to find one matches the cut
  while(out_tag==false){

    // Randomly choose an event
    rand_event = rand->Rndm() * n_event; 
    //std::cout << "Event number: " << rand_event << endl;
    data->GetMppc(rand_event,swap,newswap);

    ADC_temp.clear();
    for(int i = 0; i < fChanNum; i++) ADC_temp.push_back(data->ADC(chan_order[i]-1));

    pass_tag = EventCosmicCut(ADC_temp,ADC_cut);

    if(pass_tag==true){

      // Get the cube array
      cube_array = Get3DMatchedCubes(ADC_temp,ADC_cut);

      // 3D plot show the track (straight line)
      std::tie(track_info,chan_path,trk_graph) = Get3DTrackFit(cube_array);
      
      if(track_info[1]<40) out_tag = false;
      else out_tag = true;

    }

    if(pass_tag==true && out_tag==true){

      // Get the cube array
      //cube_array = Get3DMatchedCubes(ADC_temp,ADC_cut);

      // 3D plot show the cubes
      for(int i = 0; i < cube_array.size(); i++){
        Event_show->Fill(cube_array[i].X(),cube_array[i].Y(),cube_array[i].Z(),cube_array[i].T());
      }

      // 3D plot show the track (straight line)
      //std::tie(track_info,chan_path,trk_graph) = Get3DTrackFit(cube_array);

      // Also draw 2D MPPC on each plane
      for(int i = 0; i < fChanNum; i++){
        cube_pos = GetMPPC2DPos(chan_order[i]);
        MPPC2D_xz->Fill(cube_pos.X(),cube_pos.Z(),data->ADC(chan_order[i]-1));
        MPPC2D_yz->Fill(cube_pos.Y(),cube_pos.Z(),data->ADC(chan_order[i]-1)); 
      }
 
      //out_tag = true; 
   
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

  SetHistoSettings(Event_show);
  SetHistoSettings(bkg_temp);
  SetHistoSettings(trk_graph);

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

  TCanvas *c5 = new TCanvas("event_3d","event_3d",1200,600);
  c5->Divide(2,1);
  c5->cd(1);
  Event_show->Draw("BOX2");
  c5->cd(2);
  bkg_temp->Draw();
  trk_graph->Draw("LINE SAME");
  c5->Update();

  TString prefix = "../../../plots/scintillator_cube/";
  TString type = fPathName[file_option-1] + "/";
  TString suffix;

  suffix = prefix + type + "3Dexample_combine.png";
  c4->SaveAs(suffix);

  suffix = prefix + type + "3Dexample_combine.pdf";
  c4->SaveAs(suffix);

  suffix = prefix + type + "3Dexample_event.png";
  c5->SaveAs(suffix);

  suffix = prefix + type + "3Dexample_event.pdf";
  c5->SaveAs(suffix);

  // Save the plot into output file
  TString fout_name = "../../results/Track3DDisplay_" + fOutFileName[file_option-1] + ".root";

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();

  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c5->Write();

  Event_show->Write();
  MPPC2D_xz->Write();
  MPPC2D_yz->Write();
  trk_graph->Write();

  fout->Close();

}

// Plot the ADC distribution for each channel
void DrawMPPCLightYield(int file_option = 1, double ADC_cut = fADCCut){

  std::string fin_name;
  bool swap, newswap;

  // Input file name
  fin_name = "../../inputs/" + fFileName[file_option-1] + ".root";
  swap = fSwap[file_option-1];
  newswap = fNewSwap[file_option-1];

  TreeManager filereader(fin_name);
  Mppc *data = filereader.tmCD();

  int n_event = data->GetInputTree()->GetEntries();

  // Channel mapping choice
  if(file_option==1 || file_option==2) fChanMapChoice = 1;
  else if(file_option==3 || file_option==4) fChanMapChoice = 2;

  int chan_order[fChanNum];
  for(int i = 0; i < fChanNum; i++) chan_order[i] = fChanOrder[fChanMapChoice-1][i];

  // Create some variables
  //std::vector<double> ADC_temp;
  bool passtag;

  // Light yield distributions for each channel 0 - 17
  TString name;
  TH1D *MPPC_ly[fChanNum];
  
  for(int i = 0; i < fChanNum; i++){
    name.Form("MPPC_ly_chan%i",chan_order[i]);
    MPPC_ly[i] = new TH1D(name,name,80,0,4500);
    MPPC_ly[i]->GetXaxis()->SetTitle("ADC");
    MPPC_ly[i]->GetYaxis()->SetTitle("Number of events / bin");
    MPPC_ly[i]->SetTitle(name);
    MPPC_ly[i]->SetLineWidth(2);
    MPPC_ly[i]->SetLineColor(kBlue);
  }

  // Loop over all events
  for(int n = 0; n < n_event; n++){

    data->GetMppc(n,swap,newswap);
 
    //ADC_temp.clear();
    //for(int i = 0; i < fChanNum; i++) ADC_temp.push_back(data->ADC(chan_order[i]-1));

    //passtag = EventCosmicCut(ADC_temp,ADC_cut);
    //if(passtag==false) continue;  

    for(int i = 0; i < fChanNum; i++) MPPC_ly[i]->Fill(data->ADC(chan_order[i]-1)); 

  }

  gStyle->SetOptStat(0);

  // Drawing order array
  // MPPC channel face 1 + 3
  TCanvas *c1 = new TCanvas("MPPC2D_xz","MPPC2D_xz",1200,1200);
  c1->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c1->cd(i+1);
    MPPC_ly[i]->Draw();
    gPad->SetLogy();
  }
  c1->Update();

  // MPPC channel face 2 + 4
  TCanvas *c2 = new TCanvas("MPPC2D_yz","MPPC2D_yz",1200,1200);
  c2->Divide(3,3);
  for(int i = 0; i < 9; i++){
    c2->cd(i+1);
    MPPC_ly[i+9]->Draw();
    gPad->SetLogy();
  }
  c2->Update();

  TString prefix = "../../../plots/scintillator_cube/";
  TString type = fPathName[file_option-1] + "/";;
  TString suffix;

  suffix = prefix + type + "MPPC2D_xz.png";
  c1->SaveAs(suffix);

  suffix = prefix + type + "MPPC2D_yz.png";
  c2->SaveAs(suffix);

  // Save the plots into output file
  TString fout_name = "../../results/MPPCLightYield_" + fOutFileName[file_option-1] + ".root";

  TFile *fout = new TFile(fout_name.Data(),"recreate");
  fout->cd();
  
  for(int i = 0; i < fChanNum; i++) MPPC_ly[i]->Write();
 
  c1->Write();
  c2->Write();

  fout->Close();

} 

// ------------------------------------------------
// Auxiliary functions used in GeneralAnalysis_ZXY.cpp
// ------------------------------------------------

// Set the plot label the title
template <class T>
void SetHistoSettings(T *h){
 
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.4);
  h->GetZaxis()->SetLabelSize(0.04);
  h->GetZaxis()->SetTitleSize(0.04);
  h->GetZaxis()->SetTitleOffset(1.4);

}

// Get the channel number of nearby channels (do not consider diagonal channels)
TLorentzVector GetNearbyChannel(int i){

  TLorentzVector chan_num;

  // Old board
  if(fChanMapChoice==1){
    // Order: X = left, Y = up, Z = right, T = down
    // Face 1 + 3
    if(i==14) chan_num.SetXYZT(4,-1,-1,2);
    else if(i==4) chan_num.SetXYZT(13,-1,14,2);
    else if(i==13) chan_num.SetXYZT(-1,-1,4,3);
    else if(i==3) chan_num.SetXYZT(-1,13,12,10);
    else if(i==12) chan_num.SetXYZT(3,4,2,1);
    else if(i==2) chan_num.SetXYZT(12,14,-1,11);
    else if(i==11) chan_num.SetXYZT(1,2,-1,-1);
    else if(i==1) chan_num.SetXYZT(10,12,11,-1);
    else if(i==10) chan_num.SetXYZT(-1,3,1,-1);
    // Face 2 + 4
    else if(i==9) chan_num.SetXYZT(-1,-1,18,16);
    else if(i==18) chan_num.SetXYZT(9,-1,8,7);
    else if(i==8) chan_num.SetXYZT(18,-1,-1,17);
    else if(i==17) chan_num.SetXYZT(7,8,-1,-1);
    else if(i==7) chan_num.SetXYZT(16,18,17,15);
    else if(i==16) chan_num.SetXYZT(-1,9,7,6);
    else if(i==6) chan_num.SetXYZT(-1,16,15,-1);
    else if(i==15) chan_num.SetXYZT(6,7,5,-1);
    else if(i==5) chan_num.SetXYZT(15,17,-1,-1);
    // Not matched
    else chan_num.SetXYZT(-1,-1,-1,-1);
  }
  // New board
  else if(fChanMapChoice==2){
    // Order: X = left, Y = up, Z = right, T = down
    // Face 1 + 3
    if(i==4) chan_num.SetXYZT(-1,-1,2,17);
    else if(i==2) chan_num.SetXYZT(4,-1,6,19);
    else if(i==6) chan_num.SetXYZT(2,-1,-1,15);
    else if(i==17) chan_num.SetXYZT(-1,4,19,3);
    else if(i==19) chan_num.SetXYZT(17,2,15,1);
    else if(i==15) chan_num.SetXYZT(19,6,-1,5);
    else if(i==3) chan_num.SetXYZT(-1,17,1,-1);
    else if(i==1) chan_num.SetXYZT(3,19,5,-1);
    else if(i==5) chan_num.SetXYZT(1,15,-1,-1);
    // Face 2 + 4
    else if(i==12) chan_num.SetXYZT(-1,-1,14,23);
    else if(i==14) chan_num.SetXYZT(12,-1,8,21);
    else if(i==8) chan_num.SetXYZT(14,-1,-1,25);
    else if(i==23) chan_num.SetXYZT(-1,12,21,11);
    else if(i==21) chan_num.SetXYZT(23,14,25,13);
    else if(i==25) chan_num.SetXYZT(21,8,-1,7);
    else if(i==11) chan_num.SetXYZT(-1,23,13,-1);
    else if(i==13) chan_num.SetXYZT(11,21,7,-1);
    else if(i==7) chan_num.SetXYZT(13,25,-1,-1);
    // Not matched
    else chan_num.SetXYZT(-1,-1,-1,-1);
  }

  return chan_num;

}

// Get the two MPPC channel number corresponding to the cube
TVector3 GetCubeChannel(TVector3 cube){

  TVector3 chan_vec;

  int x = int(cube.X());
  int y = int(cube.Y());
  int z = int(cube.Z());

  // Old board
  if(fChanMapChoice==1){
  // Bottom layer
  if(x==0 && y==0 && z==0) chan_vec.SetXYZ(10,6,0);
  else if(x==1 && y==0 && z==0) chan_vec.SetXYZ(1,6,0);
  else if(x==2 && y==0 && z==0) chan_vec.SetXYZ(11,6,0);
  else if(x==0 && y==1 && z==0) chan_vec.SetXYZ(10,15,0);
  else if(x==1 && y==1 && z==0) chan_vec.SetXYZ(1,15,0);
  else if(x==2 && y==1 && z==0) chan_vec.SetXYZ(11,15,0);
  else if(x==0 && y==2 && z==0) chan_vec.SetXYZ(10,5,0);
  else if(x==1 && y==2 && z==0) chan_vec.SetXYZ(1,5,0);
  else if(x==2 && y==2 && z==0) chan_vec.SetXYZ(11,5,0);
  // Middle layer
  else if(x==0 && y==0 && z==1) chan_vec.SetXYZ(3,16,0);
  else if(x==1 && y==0 && z==1) chan_vec.SetXYZ(12,16,0); 
  else if(x==2 && y==0 && z==1) chan_vec.SetXYZ(2,16,0);
  else if(x==0 && y==1 && z==1) chan_vec.SetXYZ(3,7,0);
  else if(x==1 && y==1 && z==1) chan_vec.SetXYZ(12,7,0);
  else if(x==2 && y==1 && z==1) chan_vec.SetXYZ(2,7,0);
  else if(x==0 && y==2 && z==1) chan_vec.SetXYZ(3,17,0);
  else if(x==1 && y==2 && z==1) chan_vec.SetXYZ(12,17,0);
  else if(x==2 && y==2 && z==1) chan_vec.SetXYZ(2,17,0);
  // Top layer
  else if(x==0 && y==0 && z==2) chan_vec.SetXYZ(13,9,0);
  else if(x==1 && y==0 && z==2) chan_vec.SetXYZ(4,9,0);
  else if(x==2 && y==0 && z==2) chan_vec.SetXYZ(14,9,0);
  else if(x==0 && y==1 && z==2) chan_vec.SetXYZ(13,18,0);
  else if(x==1 && y==1 && z==2) chan_vec.SetXYZ(4,18,0);
  else if(x==2 && y==1 && z==2) chan_vec.SetXYZ(14,18,0);
  else if(x==0 && y==2 && z==2) chan_vec.SetXYZ(13,8,0);
  else if(x==1 && y==2 && z==2) chan_vec.SetXYZ(4,8,0);
  else if(x==2 && y==2 && z==2) chan_vec.SetXYZ(14,8,0);
  // other cases
  else chan_vec.SetXYZ(-1,-1,-1);
  }
  // New board
  else if(fChanMapChoice==2){
  // Bottom layer
  if(x==0 && y==0 && z==0) chan_vec.SetXYZ(3,11,0);
  else if(x==1 && y==0 && z==0) chan_vec.SetXYZ(1,11,0);
  else if(x==2 && y==0 && z==0) chan_vec.SetXYZ(5,11,0);
  else if(x==0 && y==1 && z==0) chan_vec.SetXYZ(3,13,0);
  else if(x==1 && y==1 && z==0) chan_vec.SetXYZ(1,13,0);
  else if(x==2 && y==1 && z==0) chan_vec.SetXYZ(5,13,0);
  else if(x==0 && y==2 && z==0) chan_vec.SetXYZ(3,7,0);
  else if(x==1 && y==2 && z==0) chan_vec.SetXYZ(1,7,0);
  else if(x==2 && y==2 && z==0) chan_vec.SetXYZ(5,7,0);
  // Middle layer
  else if(x==0 && y==0 && z==1) chan_vec.SetXYZ(17,23,0);
  else if(x==1 && y==0 && z==1) chan_vec.SetXYZ(19,23,0); 
  else if(x==2 && y==0 && z==1) chan_vec.SetXYZ(15,23,0);
  else if(x==0 && y==1 && z==1) chan_vec.SetXYZ(17,21,0);
  else if(x==1 && y==1 && z==1) chan_vec.SetXYZ(19,21,0);
  else if(x==2 && y==1 && z==1) chan_vec.SetXYZ(15,21,0);
  else if(x==0 && y==2 && z==1) chan_vec.SetXYZ(17,25,0);
  else if(x==1 && y==2 && z==1) chan_vec.SetXYZ(19,25,0);
  else if(x==2 && y==2 && z==1) chan_vec.SetXYZ(15,25,0);
  // Top layer
  else if(x==0 && y==0 && z==2) chan_vec.SetXYZ(4,12,0);
  else if(x==1 && y==0 && z==2) chan_vec.SetXYZ(2,12,0);
  else if(x==2 && y==0 && z==2) chan_vec.SetXYZ(6,12,0);
  else if(x==0 && y==1 && z==2) chan_vec.SetXYZ(4,14,0);
  else if(x==1 && y==1 && z==2) chan_vec.SetXYZ(2,14,0);
  else if(x==2 && y==1 && z==2) chan_vec.SetXYZ(6,14,0);
  else if(x==0 && y==2 && z==2) chan_vec.SetXYZ(4,8,0);
  else if(x==1 && y==2 && z==2) chan_vec.SetXYZ(2,8,0);
  else if(x==2 && y==2 && z==2) chan_vec.SetXYZ(6,8,0);
  // other cases
  else chan_vec.SetXYZ(-1,-1,-1);
  }

  return chan_vec;

}

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
  for(int i = 0; i <= 2; i++){
    if(ADC_temp[i]>=ADC_max){
      ADC_max = ADC_temp[i];
      index_xz = fChanOrder[fChanMapChoice-1][i];
      ADC_xz = ADC_temp[i];
    }
  }
  // Top layer (YZ plane)
  ADC_max = 0;
  for(int i = 9; i <= 11; i++){
    if(ADC_temp[i]>=ADC_max){
      ADC_max = ADC_temp[i];
      index_yz = fChanOrder[fChanMapChoice-1][i];
      ADC_yz = ADC_temp[i];
    }
  }

  cube_pos = GetMatchCubePos(index_xz,index_yz);
  cube_temp.SetXYZT(cube_pos.X(),cube_pos.Y(),cube_pos.Z(),ADC_xz+ADC_yz);
  cube_array.push_back(cube_temp);

  // Middle layer (XZ plane)
  ADC_max = 0;
  for(int i = 3; i <= 5; i++){
    if(ADC_temp[i]>=ADC_max){
      ADC_max = ADC_temp[i];
      index_xz = fChanOrder[fChanMapChoice-1][i];
      ADC_xz = ADC_temp[i];
    }
  }
  // Middle layer (YZ plane)
  ADC_max = 0;
  for(int i = 12; i <= 14; i++){
    if(ADC_temp[i]>=ADC_max){
      ADC_max = ADC_temp[i];
      index_yz = fChanOrder[fChanMapChoice-1][i];
      ADC_yz = ADC_temp[i];
    }
  }

  cube_pos = GetMatchCubePos(index_xz,index_yz);
  cube_temp.SetXYZT(cube_pos.X(),cube_pos.Y(),cube_pos.Z(),ADC_xz+ADC_yz);
  cube_array.push_back(cube_temp);

  // Bottom layer (XZ plane)
  ADC_max = 0;
  for(int i = 6; i <= 8; i++){
    if(ADC_temp[i]>=ADC_max){
      ADC_max = ADC_temp[i];
      index_xz = fChanOrder[fChanMapChoice-1][i];
      ADC_xz = ADC_temp[i];
    }
  }
  // Bottom layer (YZ plane)
  ADC_max = 0;
  for(int i = 15; i <= 17; i++){
    if(ADC_temp[i]>=ADC_max){
      ADC_max = ADC_temp[i];
      index_yz = fChanOrder[fChanMapChoice-1][i];
      ADC_yz = ADC_temp[i];
    }
  }

  cube_pos = GetMatchCubePos(index_xz,index_yz);
  cube_temp.SetXYZT(cube_pos.X(),cube_pos.Y(),cube_pos.Z(),ADC_xz+ADC_yz);
  cube_array.push_back(cube_temp);

  return cube_array;

}

// Fit the 3D track based on 3 cubes in each layer
// The output contains:
// (1) First std vector has track length, polar angle, azimuth angle and distance sum (measuring fit quality)
// (2) Second std vector has path length seen by each channel (18 in total)
// (3) A 3D track graph
std::tuple<std::vector<double>,std::vector<double>,TGraph2D*> Get3DTrackFit(std::vector<TLorentzVector> cube_array){

  // Create some parameters
  TGraph *fit_xz = new TGraph();
  TGraph *fit_yz = new TGraph();
  TF1 *fit_func = new TF1("fit_func","[0]*x+[1]");
  TF1 *fit_back;
  TGraph2D *graph = new TGraph2D();
  std::vector<double> track_info;
  std::vector<double> chan_path;
  double x_temp, y_temp;
  double ax, bx;
  double ay, by;
  double dx, dy;
  double d;
  double ang_polar, ang_azimuth;
  double dis_sum = 0;
  double a, b, c;
  double bias = fCubeSize / 2;
  double size = 3 * fCubeSize;
  double afirst, bfirst;
  double asecond, bsecond;
  TVector3 chan_pos;

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
  track_info.push_back(d); 

  // Compute the track angle with respect to positive Z axis
  ang_polar = TMath::ACos(size/d) * 180 / TMath::Pi();
  track_info.push_back(ang_polar);

  // Compute also the azimuth angle
  if(dx>=0) ang_azimuth = TMath::ATan(dy/dx) * 180 / TMath::Pi();
  else if(dx<0 && dy>0) ang_azimuth = TMath::ATan(dy/dx) * 180 / TMath::Pi() + 180;
  else ang_azimuth = TMath::ATan(dy/dx) * 180 / TMath::Pi() - 180;

  track_info.push_back(ang_azimuth);
 
  // Compute the distance sum of each projected 2D points to the fit line
  for(int i = 0; i < cube_array.size(); i++){

    // XZ plane
    a = 1;
    b = -1 * ax;
    c = -1 * bx;
    x_temp = cube_array[i].X() * fCubeSize + bias;
    y_temp = cube_array[i].Z() * fCubeSize + bias;
    dis_sum += TMath::Abs(a*x_temp+b*y_temp+c) / sqrt(a*a + b*b);
 
    // YZ plane
    a = 1;
    b = -1 * ay;
    c = -1 * by;
    x_temp = cube_array[i].Y() * fCubeSize + bias;
    y_temp = cube_array[i].Z() * fCubeSize + bias;
    dis_sum += TMath::Abs(a*x_temp+b*y_temp+c) / sqrt(a*a + b*b);

  }

  track_info.push_back(dis_sum);

  // Compute the path length seen by each channel
  double pos_xy, pos_z; // center position of the channel
  double bl[2]; // bottom left
  double br[2]; // bottom right
  double tr[2]; // top right
  double tl[2]; // top left
  double path_temp;
  double dis_temp;
  std::vector<double> intersec;
  double xy_temp, z_temp;

  for(int n = 0; n < fChanNum; n++){

    intersec.clear();

    chan_pos = GetMPPC2DPos(n+1);

    if(chan_pos.Y()==-1){ // XZ plane
      afirst = ax; bfirst = bx;
      asecond = ay; bsecond = by;
      pos_xy = chan_pos.X() * fCubeSize + bias;
      pos_z = chan_pos.Z() * fCubeSize + bias;
    }
    else{ // YZ plane
      afirst = ay; bfirst = by;
      asecond = ax; bsecond = bx;
      pos_xy = chan_pos.Y() * fCubeSize + bias;
      pos_z = chan_pos.Z() * fCubeSize + bias;
    }

    // Get positions of 4 corner points
    bl[0] = pos_xy - (fCubeSize / 2); 
    bl[1] = pos_z - (fCubeSize / 2);
    br[0] = pos_xy + (fCubeSize / 2);
    br[1] = pos_z - (fCubeSize / 2);
    tr[0] = pos_xy + (fCubeSize / 2);
    tr[1] = pos_z + (fCubeSize / 2);
    tl[0] = pos_xy - (fCubeSize / 2);
    tl[1] = pos_z + (fCubeSize / 2);

    // Now check the fit track
    // If fit track is vertical
    if(TMath::Abs(afirst)<1e-5){
      if(TMath::Abs(bfirst-pos_xy)<(fCubeSize / 2)){
        path_temp = sqrt(pow(fCubeSize,2) + pow(fCubeSize*asecond,2));
      }
      else path_temp = 0; 
    }
    // If fit track is diagonal
    else if(TMath::Abs((TMath::Abs(afirst)-1))<1e-5){
      dis_temp = TMath::Abs(pos_xy-afirst*pos_z-bfirst) / sqrt(1+afirst*afirst);
      if(dis_temp<1e-5){
        path_temp = sqrt(pow(fCubeSize,2) * 2 + pow(fCubeSize*asecond,2));
      }  
      else path_temp = 0;
    }
    // Other cases
    else{
      // Bottom line
      z_temp = bl[1];
      xy_temp = afirst * bl[1] + bfirst;
      if(xy_temp > bl[0] && xy_temp < br[0]){
        intersec.push_back(xy_temp);
        intersec.push_back(z_temp);
      }   
      // Right vertical line
      xy_temp = br[0];
      z_temp = (xy_temp - bfirst) / afirst;
      if(z_temp > br[1] && z_temp < tr[1]){
        intersec.push_back(xy_temp);
        intersec.push_back(z_temp);
      }
      // Top line
      z_temp = tr[1];
      xy_temp = afirst * tr[1] + bfirst;
      if(xy_temp > tl[0] && xy_temp < tr[0]){
        intersec.push_back(xy_temp);
        intersec.push_back(z_temp);
      }
      // Left vertical line
      xy_temp = tl[0];
      z_temp = (xy_temp - bfirst) / afirst;
      if(z_temp > bl[1] && z_temp < tl[1]){
        intersec.push_back(xy_temp);
        intersec.push_back(z_temp);
      }

      // Exactly 2 intersection points
      if(intersec.size()==4){
        path_temp = sqrt(pow(intersec[0]-intersec[2],2) + pow(intersec[1]-intersec[3],2) + pow((intersec[1]-intersec[3])*asecond,2));
      }
      else path_temp = 0;
    }

    chan_path.push_back(path_temp);

  }

  // Draw the track (straight line)
  const int npoint = 100;
  double step = size / npoint;

  for(int i = 0; i < npoint; i++){
    z_temp = step * i;
    x_temp = ax * z_temp + bx;
    y_temp = ay * z_temp + by;
    graph->SetPoint(i,x_temp,y_temp,z_temp);
  }

  return std::make_tuple(track_info,chan_path,graph);

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
  int n_midxzhit = 0;
  int n_midyzhit = 0;
  int n_tothit = 0;
  
  bool isHit; 

  // Loop over all MPPC channels (currently 18)
  for(int i = 0; i < fChanNum; i++){

    isHit = false;

    // Check if ADC higher than the cut
    if(ADC_temp[i]>=ADC_cut){
      isHit = true;
      n_tothit += 1;
    } 

    if(isHit==false) continue;

    // Old board
    if(fChanMapChoice==1){
      // Top layer (XZ plane)
      if(fChanOrder[0][i]==14 || fChanOrder[0][i]==4 || fChanOrder[0][i]==13) n_topxzhit += 1;
      // Top layer (YZ plane)
      else if(fChanOrder[0][i]==9 || fChanOrder[0][i]==18 || fChanOrder[0][i]==8) n_topyzhit += 1;
      // Bottom layer (XZ plane)
      else if(fChanOrder[0][i]==11 || fChanOrder[0][i]==1 || fChanOrder[0][i]==10) n_botxzhit += 1;
      // Bottom layer (YZ plane)
      else if(fChanOrder[0][i]==6 || fChanOrder[0][i]==15 || fChanOrder[0][i]==5) n_botyzhit += 1;
      // Middle layer (XZ plane)
      else if(fChanOrder[0][i]==3 || fChanOrder[0][i]==12 || fChanOrder[0][i]==2) n_midxzhit += 1;
      // Middle layer (YZ plane)   
      else if(fChanOrder[0][i]==17 || fChanOrder[0][i]==7 || fChanOrder[0][i]==16) n_midyzhit += 1;
    }
    // New board
    else if(fChanMapChoice==2){
      // Top layer (XZ plane)
      if(fChanOrder[1][i]==4 || fChanOrder[1][i]==2 || fChanOrder[1][i]==6) n_topxzhit += 1;
      // Top layer (YZ plane)
      else if(fChanOrder[1][i]==12 || fChanOrder[1][i]==14 || fChanOrder[1][i]==8) n_topyzhit += 1;
      // Bottom layer (XZ plane)
      else if(fChanOrder[1][i]==3 || fChanOrder[1][i]==1 || fChanOrder[1][i]==5) n_botxzhit += 1;
      // Bottom layer (YZ plane)
      else if(fChanOrder[1][i]==11 || fChanOrder[1][i]==13 || fChanOrder[1][i]==7) n_botyzhit += 1;
      // Middle layer (XZ plane)
      else if(fChanOrder[1][i]==17 || fChanOrder[1][i]==19 || fChanOrder[1][i]==15) n_midxzhit += 1;
      // Middle layer (YZ plane)   
      else if(fChanOrder[1][i]==23 || fChanOrder[1][i]==21 || fChanOrder[1][i]==25) n_midyzhit += 1;
    }
   
  }

  // Check if the requirements are satisfied 
  //if(n_tothit==6 && n_topxzhit==1 && n_topyzhit==1 && n_botxzhit==1 && n_botyzhit==1){
  if(n_topxzhit==1 && n_topyzhit==1 && n_botxzhit==1 && n_botyzhit==1){
    return true;
  }
  else{
    return false;
  }

}

// The test cubes are a 3*3*3 matrix
TVector3 GetMPPC2DPos(int MPPCChan){

  TVector3 pos;

  // Old board
  if(fChanMapChoice==1){
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
  }
  // New board
  if(fChanMapChoice==2){
  // Face 1 + 3
  if(MPPCChan==3) pos.SetXYZ(0,-1,0);
  else if(MPPCChan==1) pos.SetXYZ(1,-1,0);
  else if(MPPCChan==5) pos.SetXYZ(2,-1,0);
  else if(MPPCChan==17) pos.SetXYZ(0,-1,1);
  else if(MPPCChan==19) pos.SetXYZ(1,-1,1);
  else if(MPPCChan==15) pos.SetXYZ(2,-1,1);
  else if(MPPCChan==4) pos.SetXYZ(0,-1,2);
  else if(MPPCChan==2) pos.SetXYZ(1,-1,2);
  else if(MPPCChan==6) pos.SetXYZ(2,-1,2);
  // Face 2 + 4
  else if(MPPCChan==11) pos.SetXYZ(-1,0,0);
  else if(MPPCChan==13) pos.SetXYZ(-1,1,0);
  else if(MPPCChan==7) pos.SetXYZ(-1,2,0);
  else if(MPPCChan==23) pos.SetXYZ(-1,0,1);
  else if(MPPCChan==21) pos.SetXYZ(-1,1,1);
  else if(MPPCChan==25) pos.SetXYZ(-1,2,1);
  else if(MPPCChan==12) pos.SetXYZ(-1,0,2);
  else if(MPPCChan==14) pos.SetXYZ(-1,1,2);
  else if(MPPCChan==8) pos.SetXYZ(-1,2,2); 
  // No position available
  else pos.SetXYZ(-1,-1,-1);
  }

  return pos;

}

// Return the cube position where two fibers x and y intersect
TVector3 GetMatchCubePos(int x, int y){

  TVector3 pos;

  // Old board
  if(fChanMapChoice==1){
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
  }
  // New board
  else if(fChanMapChoice==2){
  // Bottom layer
  if(CheckIndexMatch(x,y,11,3)==true) pos.SetXYZ(0,0,0);
  else if(CheckIndexMatch(x,y,11,1)==true) pos.SetXYZ(1,0,0);
  else if(CheckIndexMatch(x,y,11,5)==true) pos.SetXYZ(2,0,0);
  else if(CheckIndexMatch(x,y,13,3)==true) pos.SetXYZ(0,1,0);
  else if(CheckIndexMatch(x,y,13,1)==true) pos.SetXYZ(1,1,0);
  else if(CheckIndexMatch(x,y,13,5)==true) pos.SetXYZ(2,1,0);
  else if(CheckIndexMatch(x,y,7,3)==true) pos.SetXYZ(0,2,0);
  else if(CheckIndexMatch(x,y,7,1)==true) pos.SetXYZ(1,2,0);
  else if(CheckIndexMatch(x,y,7,5)==true) pos.SetXYZ(2,2,0);
  // Middle layer
  else if(CheckIndexMatch(x,y,23,17)==true) pos.SetXYZ(0,0,1);
  else if(CheckIndexMatch(x,y,23,19)==true) pos.SetXYZ(1,0,1);
  else if(CheckIndexMatch(x,y,23,15)==true) pos.SetXYZ(2,0,1);
  else if(CheckIndexMatch(x,y,21,17)==true) pos.SetXYZ(0,1,1);
  else if(CheckIndexMatch(x,y,21,19)==true) pos.SetXYZ(1,1,1);
  else if(CheckIndexMatch(x,y,21,15)==true) pos.SetXYZ(2,1,1);
  else if(CheckIndexMatch(x,y,25,17)==true) pos.SetXYZ(0,2,1);
  else if(CheckIndexMatch(x,y,25,19)==true) pos.SetXYZ(1,2,1);
  else if(CheckIndexMatch(x,y,25,15)==true) pos.SetXYZ(2,2,1);
  // Top layer
  else if(CheckIndexMatch(x,y,12,4)==true) pos.SetXYZ(0,0,2);
  else if(CheckIndexMatch(x,y,12,2)==true) pos.SetXYZ(1,0,2);
  else if(CheckIndexMatch(x,y,12,6)==true) pos.SetXYZ(2,0,2);
  else if(CheckIndexMatch(x,y,14,4)==true) pos.SetXYZ(0,1,2);
  else if(CheckIndexMatch(x,y,14,2)==true) pos.SetXYZ(1,1,2);
  else if(CheckIndexMatch(x,y,14,6)==true) pos.SetXYZ(2,1,2);
  else if(CheckIndexMatch(x,y,8,4)==true) pos.SetXYZ(0,2,2);
  else if(CheckIndexMatch(x,y,8,2)==true) pos.SetXYZ(1,2,2);
  else if(CheckIndexMatch(x,y,8,6)==true) pos.SetXYZ(2,2,2);
  // Not matched
  else pos.SetXYZ(-1,-1,-1);
  }
  
  return pos;

}

bool CheckIndexMatch(int x, int y, int a, int b){

  if(x==a && y==b) return true;
  else if(x==b && y==a) return true;
  else return false;

}

// Get the index corresponding to the input channel
int ReturnIndex(int chan_num){

  int index;

  // Old board  
  if(fChanMapChoice==1){
    if(chan_num==13) index = 0;
    else if(chan_num==4) index = 1;
    else if(chan_num==14) index = 2;
    else if(chan_num==3) index = 3;
    else if(chan_num==12) index = 4;
    else if(chan_num==2) index = 5;
    else if(chan_num==10) index = 6;
    else if(chan_num==1) index = 7;
    else if(chan_num==11) index = 8;
    else if(chan_num==9) index = 9;
    else if(chan_num==18) index = 10;
    else if(chan_num==8) index = 11;
    else if(chan_num==16) index = 12;
    else if(chan_num==7) index = 13;
    else if(chan_num==17) index = 14;
    else if(chan_num==6) index = 15;
    else if(chan_num==15) index = 16;
    else if(chan_num==5) index = 17;
  }
  // New board
  else if(fChanMapChoice==2){
    if(chan_num==4) index = 0;
    else if(chan_num==2) index = 1;
    else if(chan_num==6) index = 2;
    else if(chan_num==17) index = 3;
    else if(chan_num==19) index = 4;
    else if(chan_num==15) index = 5;
    else if(chan_num==3) index = 6;
    else if(chan_num==1) index = 7;
    else if(chan_num==5) index = 8;
    else if(chan_num==12) index = 9;
    else if(chan_num==14) index = 10;
    else if(chan_num==8) index = 11;
    else if(chan_num==23) index = 12;
    else if(chan_num==21) index = 13;
    else if(chan_num==25) index = 14;
    else if(chan_num==11) index = 15;
    else if(chan_num==13) index = 16;
    else if(chan_num==7) index = 17;
  }

  return index;

}

bool CheckTrackVertical(std::vector<TLorentzVector> cube_array){

  bool x_match = true;
  bool y_match = true;
  
  for(int n = 1; n < cube_array.size(); n++){
    
    if(cube_array[n].X() != cube_array[n-1].X()) x_match = false;
    
    if(cube_array[n].Y() != cube_array[n-1].Y()) y_match = false;
    
  }
  
  bool ver_tag;
  
  if(x_match == true && y_match == true) ver_tag = true;
  else ver_tag = false;
  
  return ver_tag;

}
