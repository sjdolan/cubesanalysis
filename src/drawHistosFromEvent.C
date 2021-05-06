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

#include "EventAnalysis.cpp"

bool isSFGDChannelSwap=false;

//TreeManager filereader("/eos/home-s/sdolan/cubesWork/analysis/fromUmut/GluedCubes/SpecialRun_GluedCubesNoTeflon_15April2021.root");
//TreeManager filereader("/eos/home-s/sdolan/cubesWork/analysis/fromUmut/GluedCubes/SpecialRun_GluedCubesNoTeflon_14April2021.root");

//TreeManager filereader("/eos/home-d/dsgalabe/3DprintScint/RD/Characterization_April_2021/Data/CosmicMuon_GluedCubesandFibersWithTeflonLayers_01May2021.root");

// Combined data with Teflon before gluing the fibres into the connectors
//TreeManager filereader("/eos/home-s/sdolan/cubesWork/data/SpecialRun_GluedCubesWithTeflon_19to21AprilHadded.root");

// New data after the gluing. Around 48 hours of data taking
//TreeManager filereader("/eos/home-s/sdolan/cubesWork/data/CosmicMuon_GluedCubesandFibersWithTeflonLayers_01To02May_hadded.root");
//string outPlotDir = "/eos/home-s/sdolan/cubesWork/git/plots/Glue_";
// SFGD cubes
TreeManager filereader("/eos/home-s/sdolan/cubesWork/data/SpecialRun_SuperFGDCubes_04to06MayHadded.root");
string outPlotDir = "/eos/home-s/sdolan/cubesWork/git/plots/SFGD_";

using namespace std;

//gStyle->SetPaintTextFormat("6.1f");

//Function prototypes
std::vector<TH2D*> twoFaces(Long64_t, double, double, bool);
TH3D* threeDimWLSF(Long64_t, double, double, bool);
TH3D* threeDimMPPC(Long64_t, double, double, bool);
TH1D* totalLightYield();
TVector3 getMPPCPos(int);
void showAll(Long64_t, double, double);
std::vector<double> getPathLengthThroughChannel(Long64_t, UInt_t, bool);
double getPolarAngle(Long64_t, UInt_t, bool);
bool passCosmicSelection(Long64_t, UInt_t, bool);
double GetMPVPos(TH1D*);


void setSFGDChannelSwap(){
  isSFGDChannelSwap=true; // Channels 8 and 18 are swapped for SFGD data
}



double GetMPVPos(TH1D* h){
    int binmax = h->GetMaximumBin();
    double x = h->GetXaxis()->GetBinCenter(binmax);
    return x;
}


void showHighLYEvent(UInt_t lyLim, int seed=0){
    //Initalise RNG generator
    TRandom3* rando = new TRandom3();
    rando->SetSeed(seed);
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    int randEvNum = 0;

    UInt_t totalLY = 0;
    while(true){
        randEvNum = rando->Rndm()*entries;
        data->GetMppc(randEvNum, isSFGDChannelSwap);
        totalLY = 0;
        for (int i=0;i<18;i++) totalLY += data->ADC(i);
        if(totalLY>lyLim) break;
    }

    std::cout << "Chosen event number: " << randEvNum << std::endl;
    std::cout << "  Light Yield (ADC): " << totalLY << std::endl;

    showAll(randEvNum, 1, 1000);
}

void showHighADCCoincidenceEvent(UInt_t lyLim=1000, UInt_t hitLim=3, int seed=0){
    //Initalise RNG generator
    TRandom3* rando = new TRandom3();
    rando->SetSeed(seed);
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    int randEvNum = 0;

    UInt_t totalLY = 0;
    UInt_t nHits = 0;
    while(true){
        randEvNum = rando->Rndm()*entries;
        data->GetMppc(randEvNum, isSFGDChannelSwap);
        totalLY = 0;
        nHits = 0;
        for (int i=0;i<18;i++){
            totalLY += data->ADC(i);
            if(data->ADC(i) > lyLim) nHits+=1;
        }
        if(nHits>hitLim) break;
    }

    std::cout << "Chosen event number: " << randEvNum << std::endl;
    std::cout << "  Light Yield (ADC): " << totalLY << std::endl;
    std::cout << "     Number of hits: " << nHits << std::endl;

    showAll(randEvNum, 1, 1000);
}

void showCosmicCoincidenceEvent(UInt_t lyLim=1000, int seed=0){
    //Initalise RNG generator
    TRandom3* rando = new TRandom3();
    rando->SetSeed(seed);
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    Long64_t randEvNum = 0;

    UInt_t totalLY = 0;
    UInt_t nHits = 0;
    UInt_t nHitsTopX = 0;
    UInt_t nHitsTopY = 0;
    UInt_t nHitsBotX = 0;
    UInt_t nHitsBotY = 0;
    bool isHit = false;
    std::vector<UInt_t> hitMPPC; 
    std::vector<UInt_t> hitADC; 

    while(true){
        randEvNum = rando->Rndm()*entries;
        if( passCosmicSelection(randEvNum, lyLim, true) ) break;
    }
    showAll(randEvNum, 1, 1000);
}


void showAll(Long64_t evNum, double limLow=1, double limHigh=300){

    std::vector<TH2D*> listOfHist = twoFaces(evNum, limLow, limHigh, true);
    TH3D* histoWLSF = threeDimWLSF(evNum, limLow, limHigh, true);
    TH3D* histoMPPC = threeDimMPPC(evNum, limLow, limHigh, true);

    gStyle->SetOptTitle(1); 
    TCanvas *c1 = new TCanvas("c1","multipads",1350,950);
    gStyle->SetOptStat(0);
    c1->Divide(2,2);
    c1->cd(1);
    listOfHist[0]->Draw("colzTEXT45");
    c1->cd(2);
    listOfHist[1]->Draw("colzTEXT45");
    c1->cd(3);
    histoWLSF->Draw("BOX2");
    c1->cd(4);
    histoMPPC->Draw("BOX2");

}



std::vector<TH2D*> twoFaces(Long64_t evNum, double limLow=1, double limHigh=300, bool noDraw=false){

    //Setup hhistos
    TH2D* facexzHisto = new TH2D("xz", "XZ View;X [cm]; Z [cm]", 7,0,7, 7,0,7);
    TH2D* faceyzHisto = new TH2D("yz", "YZ View;Y [cm]; Z [cm]", 7,0,7, 7,0,7);

    facexzHisto->GetZaxis()->SetRangeUser(limLow, limHigh);
    faceyzHisto->GetZaxis()->SetRangeUser(limLow, limHigh);

	Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    data->GetMppc(evNum, isSFGDChannelSwap);
    std::cout << "Loading event number: " << evNum << std::endl;

    //Loop over the MPPCs in the event
    for(int i=1; i<=18; i++){
        TVector3 pos = getMPPCPos(i);
        //Face 1 or Face 3
        if(pos.Y()==0 || pos.Y()==7){
            facexzHisto->Fill(pos.X(), pos.Z(), (double)data->ADC(i-1));
        }
        //Face 2 or Face 4
        else if(pos.X()==0 || pos.X()==7){
            faceyzHisto->Fill(pos.Y(), pos.Z(), (double)data->ADC(i-1));
        }
        else{
            std::cout << "WARNING: found invalid MPPC position for MPPC " << i << std::endl;
            pos.Print();
        }
    }

    if(!noDraw){
        TCanvas *c1 = new TCanvas("c1","multipads",1600,700);
        gStyle->SetOptStat(0);
        c1->Divide(2,1);
        c1->cd(1);
        facexzHisto->Draw("colzTEXT45");
        c1->cd(2);
        faceyzHisto->Draw("colzTEXT45");
    }

    std::vector<TH2D*> listOfHist; 
    listOfHist.push_back(facexzHisto);
    listOfHist.push_back(faceyzHisto);

    return(listOfHist);

}

TH3D* threeDimWLSF(Long64_t evNum, double limLow=1, double limHigh=300, bool noDraw=false){
    //Setup hhistos
    TH3D* histo = new TH3D("xyzWLSF", "3D WLSF View;X [mm]; Y [mm]; Z[mm]", 14,0,70, 14,0,70 , 14,0,70);

    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    data->GetMppc(evNum, isSFGDChannelSwap);
    std::cout << "Loading event number: " << evNum << std::endl;

    //Loop over the MPPCs in the event
    for(int i=1; i<=18; i++){
        TVector3 pos = getMPPCPos(i);
        //Face 1 or Face 3
        if(pos.Y()==0 || pos.Y()==7){
            for(int j=2; j<=12; j++){
              histo->Fill(10*pos.X()+5, 5*j, 10*pos.Z()+5, (double)data->ADC(i-1));
            }
        }
        //Face 2 or Face 4
        else if(pos.X()==0 || pos.X()==7){
            for(int j=2; j<=12; j++){
              histo->Fill(5*j, 10*pos.Y()+5, 10*pos.Z()+5, (double)data->ADC(i-1));
            }
        }
        else{
            std::cout << "WARNING: found invalid MPPC position for MPPC " << i << std::endl;
            pos.Print();
        }
    }

    if(!noDraw) histo->Draw("BOX2");

    return histo;
}


TH3D* threeDimMPPC(Long64_t evNum, double limLow=1, double limHigh=300, bool noDraw=false){
    //Setup hhistos
    TH3D* histo = new TH3D("xyzMPPC", "3D MPPC View;X [cm]; Y [cm]; Z[cm]", 7,0,7, 7,0,7 , 7,0,7);

    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    data->GetMppc(evNum, isSFGDChannelSwap);
    std::cout << "Loading event number: " << evNum << std::endl;

    //Loop over the MPPCs in the event
    for(int i=1; i<=18; i++){
        TVector3 pos = getMPPCPos(i);
        //Face 1 or Face 3
        if(pos.Y()==0 || pos.Y()==7){
            histo->Fill(pos.X(), 0.0, pos.Z(), (double)data->ADC(i-1));
        }
        //Face 2 or Face 4
        else if(pos.X()==0 || pos.X()==7){
            histo->Fill(0.0, pos.Y(), pos.Z(), (double)data->ADC(i-1));
        }
        else{
            std::cout << "WARNING: found invalid MPPC position for MPPC " << i << std::endl;
            pos.Print();
        }
    }

    if(!noDraw) histo->Draw("BOX2");

    return histo;
}

TH1D* totalLightYield()
{
	Mppc *data = filereader.tmCD(); 
	Long64_t entries = data->GetInputTree()->GetEntries();

    TH1D* lyHisto = new TH1D("totalLightYield", "totalLightYield", 100, 0, 4000);
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry, isSFGDChannelSwap);
        double totalLY = 0;
        for (int i=0;i<18;i++) totalLY += data->ADC(i);
        lyHisto->Fill(totalLY);
	}
    TCanvas* c1 = new TCanvas();
    lyHisto->Draw();

    return lyHisto;
}

std::vector<TH1D*> perMPPCLightYield()
{
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();

    std::vector<TH1D*> listOfHist; 
    for (int i=1;i<=18;i++){
        TH1D* lyHisto = new TH1D(Form("lyHisto%d",i), Form("MPPC %d; Light Yield [ADC]; Counts",i), 300, 0, 4000);
        listOfHist.push_back(lyHisto);
    }

    for( Long64_t ientry=0; ientry<entries; ++ientry )
    {
        data->GetMppc(ientry, isSFGDChannelSwap);
        for (int i=0;i<18;i++){
            listOfHist[i]->Fill(data->ADC(i));
        }
    }
    gStyle->SetOptTitle(1); 
    TCanvas *c1 = new TCanvas("c1","c1",1900,950);
    gStyle->SetOptStat(0);
    c1->Divide(6,3);
    for (int i=0;i<18;i++){
        gPad->SetLogy();
        c1->cd(i+1);
        listOfHist[i]->Draw();
    }

    return listOfHist;
}

std::vector<TH1D*> perMPPCCosmicLightYieldAndHits(UInt_t lyLim=1500, bool applyPathLengthCorr=false, bool onlyVert=false, int nCosmicsLimit=0)
{
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();


    //Setup histos
    TH1D* hitsHisto = new TH1D("hitsHisto", "hitsHisto;MPPC Number; Hits From Cosmics", 18, 1, 19);

    TH1D* LYDevHisto = new TH1D("LYDevHisto", "LYDevHisto;MPPC Number; Deviation From Average LY [%]", 18, 1, 19);
    TH1D* meanLYHisto = new TH1D("meanLYHisto", "meanLYHisto;MPPC Number; Average LY", 18, 1, 19);

    TH1D* LYDevHisto13 = new TH1D("LYDevHisto13", "Faces 1 and 3;Offset MPPC Number; Deviation From Average LY [%]", 9, 1, 10);
    TH1D* meanLYHisto13 = new TH1D("meanLYHisto13", "Faces 1 and 3;Offset MPPC Number; Average LY", 9, 1, 10);

    TH1D* LYDevHisto24 = new TH1D("LYDevHisto24", "Faces 2 and 4;Offset MPPC Number; Deviation From Average LY [%]", 9, 1, 10);
    TH1D* meanLYHisto24 = new TH1D("meanLYHisto24", "Faces 2 and 4;Offset MPPC Number; Average LY", 9, 1, 10);

    TH2D* facexzHisto = new TH2D("MeanCosmicLY_XZ", "Mean Cosmic LY;X [cm]; Z [cm]", 7,0,7, 7,0,7);
    TH2D* faceyzHisto = new TH2D("MeanCosmicLY_YZ", "Mean Cosmic LY;Y [cm]; Z [cm]", 7,0,7, 7,0,7);

    TH1D* MPVLYDevHisto = new TH1D("MPVLYDevHisto", "MPVLYDevHisto;MPPC Number; Deviation From Average LY [%]", 18, 1, 19);
    TH1D* MPVLYHisto = new TH1D("MPVLYHisto", "MPVLYHisto;MPPC Number; Average LY", 18, 1, 19);

    TH1D* MPVLYDevHisto13 = new TH1D("MPVLYDevHisto13", "Faces 1 and 3;Offset MPPC Number; Deviation From Average LY [%]", 9, 1, 10);
    TH1D* MPVLYHisto13 = new TH1D("MPVLYHisto13", "Faces 1 and 3;Offset MPPC Number; Average LY", 9, 1, 10);

    TH1D* MPVLYDevHisto24 = new TH1D("MPVLYDevHisto24", "Faces 2 and 4;Offset MPPC Number; Deviation From Average LY [%]", 9, 1, 10);
    TH1D* MPVLYHisto24 = new TH1D("MPVLYHisto24", "Faces 2 and 4;Offset MPPC Number; Average LY", 9, 1, 10);

    TH2D* facexzMPVHisto = new TH2D("MPVCosmicLY_XZ", "MPV Cosmic LY;X [cm]; Z [cm]", 7,0,7, 7,0,7);
    TH2D* faceyzMPVHisto = new TH2D("MPVCosmicLY_YZ", "MPV Cosmic LY;Y [cm]; Z [cm]", 7,0,7, 7,0,7);

    facexzHisto->GetZaxis()->SetRangeUser(lyLim, 4000);
    faceyzHisto->GetZaxis()->SetRangeUser(lyLim, 4000);

    facexzMPVHisto->GetZaxis()->SetRangeUser(lyLim, 4000);
    faceyzMPVHisto->GetZaxis()->SetRangeUser(lyLim, 4000);

    std::vector<double> pathLength;
    std::vector<TH1D*> listOfHist; 
    for (int i=1;i<=18;i++){
        TH1D* lyHisto = new TH1D(Form("lyHisto%d",i), Form("MPPC %d; Light Yield [ADC]; Counts",i), 30, lyLim, 4000);
        listOfHist.push_back(lyHisto);
        pathLength.push_back(10.0);
    }

    int nCosmicsFound = 0;

    for( Long64_t ientry=0; ientry<entries; ++ientry )
    {
        // Check if we pass the cosmic selection for this event
        if( !passCosmicSelection(ientry, lyLim, false) ) continue;

        if(onlyVert){
            double polarAng = getPolarAngle(ientry, lyLim, true);
            if (abs(polarAng)>1E-5) continue;
        }

        nCosmicsFound++;
        if(nCosmicsLimit!=0 && nCosmicsFound>nCosmicsLimit) break;

        if(applyPathLengthCorr){
            std::cout << "Processing event " << ientry << " out of " << entries << std::endl;
            pathLength.clear();
            pathLength = getPathLengthThroughChannel(ientry, lyLim, true);
        }

        data->GetMppc(ientry, isSFGDChannelSwap);
        for (int i=0;i<18;i++){
            if(data->ADC(i) > lyLim){
                hitsHisto->Fill(i+1);
                // Only apply the correction if the MPPC did not saturate
                if(data->ADC(i) < 4000.) listOfHist[i]->Fill(data->ADC(i)*(10.0/pathLength[i]));
                else                     listOfHist[i]->Fill(data->ADC(i));
                //std::cout << "Path lenth for MPPC " << i << " is " << pathLength[i] << std::endl;
            }
        }
    }



    //Loop over the MPPCs 
    for(int i=1; i<=18; i++){
        TVector3 pos = getMPPCPos(i);
        //Face 1 or Face 3
        if(pos.Y()==0 || pos.Y()==7){
            facexzHisto->Fill(pos.X(), pos.Z(), listOfHist[i-1]->GetMean());
            facexzMPVHisto->Fill(pos.X(), pos.Z(), GetMPVPos(listOfHist[i-1]));
        }
        //Face 2 or Face 4
        else if(pos.X()==0 || pos.X()==7){
            faceyzHisto->Fill(pos.Y(), pos.Z(), listOfHist[i-1]->GetMean());
            faceyzMPVHisto->Fill(pos.Y(), pos.Z(), GetMPVPos(listOfHist[i-1]));

        }
        else{
            std::cout << "WARNING: found invalid MPPC position for MPPC " << i << std::endl;
            pos.Print();
        }
    }

    //Find the average light yield across all MPPCs
    double meanOfMeans = 0;
    double meanOfMeans13 = 0;
    double meanOfMeans24 = 0;
    for(int i=1; i<=18; i++){
        meanOfMeans += listOfHist[i-1]->GetMean();
        if(i<5 || (i>9 && i<15)) meanOfMeans13 += listOfHist[i-1]->GetMean(); // Face 1 or 3
        else meanOfMeans24 += listOfHist[i-1]->GetMean(); // Face 2 or 4

    }
    meanOfMeans = meanOfMeans / 18.0;
    meanOfMeans13 = meanOfMeans13 / 9.0;
    meanOfMeans24 = meanOfMeans24 / 9.0;

    //Find the average of the MPV value of the light yield across all MPPCs
    double meanOfMPVs = 0;
    double meanOfMPVs13 = 0;
    double meanOfMPVs24 = 0;
    for(int i=1; i<=18; i++){
        meanOfMPVs += GetMPVPos(listOfHist[i-1]);
        if(i<5 || (i>9 && i<15)) meanOfMPVs13 += GetMPVPos(listOfHist[i-1]); // Face 1 or 3
        else meanOfMPVs24 += GetMPVPos(listOfHist[i-1]); // Face 2 or 4

    }
    meanOfMPVs = meanOfMPVs / 18.0;
    meanOfMPVs13 = meanOfMPVs13 / 9.0;
    meanOfMPVs24 = meanOfMPVs24 / 9.0;

    int face13Count=1;
    int face24Count=1;
    for(int i=1; i<=18; i++){
        double avLY = listOfHist[i-1]->GetMean();
        double MPVLY = GetMPVPos(listOfHist[i-1]);

        meanLYHisto->SetBinContent(i, avLY);
        LYDevHisto->SetBinContent(i, 100.0*(avLY-meanOfMeans)/meanOfMeans);
        MPVLYHisto->SetBinContent(i, MPVLY);
        MPVLYDevHisto->SetBinContent(i, 100.0*(MPVLY-meanOfMPVs)/meanOfMPVs);

        //Face 1 or Face 3
        if(i<5 || (i>9 && i<15)){
            meanLYHisto13->SetBinContent(face13Count, avLY);
            LYDevHisto13->SetBinContent(face13Count, 100.0*(avLY-meanOfMeans13)/meanOfMeans13);
            MPVLYHisto13->SetBinContent(face13Count, MPVLY);
            MPVLYDevHisto13->SetBinContent(face13Count, 100.0*(MPVLY-meanOfMPVs13)/meanOfMPVs13);
            face13Count++;
        }
        //Face 2 or Face 4
        else{
            meanLYHisto24->SetBinContent(face24Count, avLY);
            LYDevHisto24->SetBinContent(face24Count, 100.0*(avLY-meanOfMeans24)/meanOfMeans24);
            MPVLYHisto24->SetBinContent(face24Count, MPVLY);
            MPVLYDevHisto24->SetBinContent(face24Count, 100.0*(MPVLY-meanOfMPVs24)/meanOfMPVs24);
            face24Count++;
        }
    }

    std::cout << "Average of average light yields is " << meanOfMeans << " [ADC]" << std::endl;
    std::cout << "Average of average light yields on faces 1 and 3 is " << meanOfMeans13 << " [ADC]" << std::endl;
    std::cout << "Average of average light yields on faces 2 and 4 is " << meanOfMeans24 << " [ADC]" << std::endl;

    std::cout << "Average of MPV light yields is " << meanOfMPVs << " [ADC]" << std::endl;
    std::cout << "Average of MPV light yields on faces 1 and 3 is " << meanOfMPVs13 << " [ADC]" << std::endl;
    std::cout << "Average of MPV light yields on faces 2 and 4 is " << meanOfMPVs24 << " [ADC]" << std::endl;

    gStyle->SetOptTitle(1); 
    TCanvas *c1 = new TCanvas("c1","c1",1900,950);
    gStyle->SetOptStat(0);
    c1->Divide(6,3);
    for (int i=0;i<18;i++){
        c1->cd(i+1);
        listOfHist[i]->Draw();
    }
    c1->SaveAs(Form("%sLYHistos.png", outPlotDir.c_str()));


    TCanvas *c2 = new TCanvas();
    hitsHisto->Draw();
    c2->SaveAs(Form("%shitsHisto.png", outPlotDir.c_str()));


    TCanvas *c3 = new TCanvas("multipads","multipads",1600,700);
    gStyle->SetOptStat(0);
    c3->Divide(2,1);
    c3->cd(1);
    facexzHisto->Draw("colzTEXT45");
    c3->cd(2);
    faceyzHisto->Draw("colzTEXT45");
    c3->SaveAs(Form("%smeanLY2DHistos.png", outPlotDir.c_str()));


    TCanvas *c4 = new TCanvas();
    meanLYHisto->Draw();
    c4->SaveAs(Form("%smeanLYHisto.png", outPlotDir.c_str()));


    TCanvas *c5 = new TCanvas();
    LYDevHisto->GetYaxis()->SetRangeUser(-50.0, 50.0);
    LYDevHisto->Draw();
    c5->SaveAs(Form("%sLYDevHisto.png", outPlotDir.c_str()));


    TCanvas *c6 = new TCanvas("face13","face13",1600,700);
    gStyle->SetOptStat(0);
    c6->Divide(2,1);
    c6->cd(1);
    meanLYHisto13->Draw();
    c6->cd(2);
    LYDevHisto13->GetYaxis()->SetRangeUser(-50.0, 50.0);
    LYDevHisto13->Draw();
    c6->SaveAs(Form("%sface13Histos.png", outPlotDir.c_str()));


    TCanvas *c7 = new TCanvas("face24","face24",1600,700);
    gStyle->SetOptStat(0);
    c7->Divide(2,1);
    c7->cd(1);
    meanLYHisto24->Draw();
    c7->cd(2);
    LYDevHisto24->GetYaxis()->SetRangeUser(-50.0, 50.0);
    LYDevHisto24->Draw();
    c7->SaveAs(Form("%sface24Histos.png", outPlotDir.c_str()));


    // MPV Canvas

    TCanvas *c3MPV = new TCanvas("multipadsMPV","multipadsMPV",1600,700);
    gStyle->SetOptStat(0);
    c3MPV->Divide(2,1);
    c3MPV->cd(1);
    facexzMPVHisto->Draw("colzTEXT45");
    c3MPV->cd(2);
    faceyzMPVHisto->Draw("colzTEXT45");
    c3MPV->SaveAs(Form("%sMPVLY2DHistos.png", outPlotDir.c_str()));

    TCanvas *c4MPV = new TCanvas();
    MPVLYHisto->Draw();
    c4MPV->SaveAs(Form("%sMPVLYHisto.png", outPlotDir.c_str()));

    TCanvas *c5MPV = new TCanvas();
    MPVLYDevHisto->GetYaxis()->SetRangeUser(-50.0, 50.0);
    MPVLYDevHisto->Draw();
    c5MPV->SaveAs(Form("%sMPVLYDevHisto.png", outPlotDir.c_str()));

    TCanvas *c6MPV = new TCanvas("face13MPV","face13MPV",1600,700);
    gStyle->SetOptStat(0);
    c6MPV->Divide(2,1);
    c6MPV->cd(1);
    MPVLYHisto13->Draw();
    c6MPV->cd(2);
    MPVLYDevHisto13->GetYaxis()->SetRangeUser(-50.0, 50.0);
    MPVLYDevHisto13->Draw();
    c6MPV->SaveAs(Form("%sMPVface13Histos.png", outPlotDir.c_str()));


    TCanvas *c7MPV = new TCanvas("face24MPV","face24MPV",1600,700);
    gStyle->SetOptStat(0);
    c7MPV->Divide(2,1);
    c7MPV->cd(1);
    MPVLYHisto24->Draw();
    c7MPV->cd(2);
    MPVLYDevHisto24->GetYaxis()->SetRangeUser(-50.0, 50.0);
    MPVLYDevHisto24->Draw();;
    c7MPV->SaveAs(Form("%sMPVface24Histos.png", outPlotDir.c_str()));


    return listOfHist;
}



TH1D* nHitsHisto(UInt_t lyLim=1500)
{
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();

    TH1D* nHitsHisto = new TH1D("nHitsHisto", "nHitsHisto", 9, 0, 9);
    UInt_t nHits = 0;
    for( Long64_t ientry=0; ientry<entries; ++ientry )
    {
        data->GetMppc(ientry, isSFGDChannelSwap);
        UInt_t totalLY = 0;
        nHits = 0;
        for (int i=0;i<18;i++){
            totalLY += data->ADC(i);
            if(data->ADC(i) > lyLim) nHits+=1;
        }
        nHitsHisto->Fill((double)nHits);
    }
    TCanvas* c1 = new TCanvas();
    nHitsHisto->Draw();

    return nHitsHisto;
}

TH1D* coincidenceRateHisto(UInt_t lyLim=1500)
{
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();

    TH1D* isCoincidence = new TH1D("isCoincidence", "isCoincidence", 2, 0, 2);

    for( Long64_t ientry=0; ientry<entries; ++ientry )
    {
        if(passCosmicSelection(ientry, lyLim, false)) isCoincidence->Fill(1.5);
        else isCoincidence->Fill(0.5);
    }
    TCanvas* c1 = new TCanvas();
    isCoincidence->Draw();
    isCoincidence->Print("all");

    return isCoincidence;
}

std::vector<double> getPathLengthThroughChannel(Long64_t evNum, UInt_t lyLim=1500, bool noDraw=false){

    std::vector<double> ADC_temp;
    std::vector<TLorentzVector> cube_array;

    TGraph2D *trk_graph;
    std::vector<double> track_info;
    std::vector<double> chan_path;

    Mppc *data = filereader.tmCD();
    data->GetMppc(evNum, isSFGDChannelSwap);
    for(int i = 0; i < 18; i++) ADC_temp.push_back(data->ADC(i));

    // Get the positions (+ light yields) of matched cubes
    // Each layer only consider the highest ADC in each plane, and use the center position of cube
    cube_array = Get3DMatchedCubes(ADC_temp,lyLim);

    // Get the total length of the track and angle with respect to positive Z axis
    std::tie(track_info,chan_path,trk_graph) = Get3DTrackFit(cube_array);

    if(!noDraw){
        double trk_len = track_info[0];
        double trk_ang_pol = track_info[1];
        double trk_ang_azi = track_info[2];

        TH1D* pathLengthHisto = new TH1D("pathLengthHisto", "pathLengthHisto;MPPC Number; Path length in event", 18, 1, 19);
        for(int i = 0; i < 18; i++) pathLengthHisto->SetBinContent(i+1, chan_path[i]);
        
        std::cout << "Track lenth is: " << trk_len << std::endl;
        std::cout << "Track polar angle is: " << trk_ang_pol << std::endl;
        std::cout << "Track azimuthal angle is: " << trk_ang_azi << std::endl;

        TCanvas* c1 = new TCanvas();
        pathLengthHisto->Draw();

    }

    return chan_path;
}

double getPolarAngle(Long64_t evNum, UInt_t lyLim=1500, bool noDraw=false){

    std::vector<double> ADC_temp;
    std::vector<TLorentzVector> cube_array;

    TGraph2D *trk_graph;
    std::vector<double> track_info;
    std::vector<double> chan_path;

    Mppc *data = filereader.tmCD();
    data->GetMppc(evNum, isSFGDChannelSwap);
    for(int i = 0; i < 18; i++) ADC_temp.push_back(data->ADC(i));

    // Get the positions (+ light yields) of matched cubes
    // Each layer only consider the highest ADC in each plane, and use the center position of cube
    cube_array = Get3DMatchedCubes(ADC_temp,lyLim);

    // Get the total length of the track and angle with respect to positive Z axis
    std::tie(track_info,chan_path,trk_graph) = Get3DTrackFit(cube_array);

    double trk_len = track_info[0];
    double trk_ang_pol = track_info[1];
    double trk_ang_azi = track_info[2];

    return trk_ang_pol;
}



bool passCosmicSelection(Long64_t evNum, UInt_t lyLim, bool verbose){
    bool pass = false;

    Mppc *data = filereader.tmCD(); 
    data->GetMppc(evNum, isSFGDChannelSwap);

    UInt_t totalLY = 0;
    UInt_t nHits = 0;
    UInt_t nHitsTopX = 0;
    UInt_t nHitsTopY = 0;
    UInt_t nHitsBotX = 0;
    UInt_t nHitsBotY = 0;
    bool isHit = false;
    std::vector<UInt_t> hitMPPC; 
    std::vector<UInt_t> hitADC; 

    for (Int_t i=0;i<18;i++){
        totalLY += data->ADC(i);
        if(data->ADC(i) > lyLim) isHit=true;
        else isHit = false;
        Int_t mppcNum = i+1;
        if(isHit){
            if( mppcNum == 4 || mppcNum == 14 || mppcNum == 13) nHitsTopX+=1;
            if( mppcNum == 9 || mppcNum == 8  || mppcNum == 18) nHitsTopY+=1;
            if( mppcNum == 1 || mppcNum == 11 || mppcNum == 10) nHitsBotX+=1;
            if( mppcNum == 6 || mppcNum == 5  || mppcNum == 15) nHitsBotY+=1;
            nHits++;
            hitMPPC.push_back(mppcNum);
            hitADC.push_back(data->ADC(i));
        }
    }
    if(nHits==6 && nHitsTopX==1 && nHitsTopY==1 && nHitsBotX==1 && nHitsBotY==1) pass = true;
    //if(nHitsTopX==1 && nHitsTopY==1 && nHitsBotX==1 && nHitsBotY==1) pass = true;
    else pass = false;

    if(pass == true && verbose){
        std::cout << "Chosen event number: " << evNum << std::endl;
        std::cout << "  Light Yield (ADC): " << totalLY << std::endl;
        std::cout << "     Number of hits: " << nHits << std::endl;
        for (int i=0; i<hitMPPC.size(); i++){
          std::cout << "Hit MPPC: " << hitMPPC[i] << ", ADC counts: " << hitADC[i] << std::endl;
        }
    }

    return pass;
}

TVector3 getMPPCPos(int mppcNum){
    TVector3 pos;
    double offset = 0.5; 

    //Face 1
    if(mppcNum == 1)  pos.SetXYZ(3,0,1); 
    else if(mppcNum == 2)  pos.SetXYZ(4,0,2); 
    else if(mppcNum == 3)  pos.SetXYZ(2,0,2); 
    else if(mppcNum == 4)  pos.SetXYZ(3,0,3); 
    //Face 2
    else if(mppcNum == 5)  pos.SetXYZ(0,2,1+offset); 
    else if(mppcNum == 6)  pos.SetXYZ(0,4,1+offset); 
    else if(mppcNum == 7)  pos.SetXYZ(0,3,2+offset); 
    else if(mppcNum == 8)  pos.SetXYZ(0,2,3+offset);
    else if(mppcNum == 9)  pos.SetXYZ(0,4,3+offset); 
    //Face 3 
    else if(mppcNum == 10) pos.SetXYZ(2,7,1); 
    else if(mppcNum == 11) pos.SetXYZ(4,7,1); 
    else if(mppcNum == 12) pos.SetXYZ(3,7,2); 
    else if(mppcNum == 13) pos.SetXYZ(2,7,3); 
    else if(mppcNum == 14) pos.SetXYZ(4,7,3); 
    //Face 4
    else if(mppcNum == 15) pos.SetXYZ(7,3,1+offset); 
    else if(mppcNum == 16) pos.SetXYZ(7,4,2+offset); 
    else if(mppcNum == 17) pos.SetXYZ(7,2,2+offset); 
    else if(mppcNum == 18) pos.SetXYZ(7,3,3+offset); 
    else pos.SetXYZ(-1,-1,-1);

    return pos;
}
