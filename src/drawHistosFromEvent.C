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

//TreeManager filereader("/eos/home-s/sdolan/cubesWork/analysis/fromUmut/GluedCubes/SpecialRun_GluedCubesNoTeflon_15April2021.root");
//TreeManager filereader("/eos/home-s/sdolan/cubesWork/analysis/fromUmut/GluedCubes/SpecialRun_GluedCubesNoTeflon_14April2021.root");
TreeManager filereader("/eos/user/d/dsgalabe/3DprintScint/RD/Characterization_April_2021/Data/GluedCubes/SpecialRun_GluedCubesWithTeflon_19April2021.root");

using namespace std;

//gStyle->SetPaintTextFormat("6.1f");

//Function prototypes
vector<TH2D*> twoFaces(Long64_t, double, double, bool);
TH3D* threeDimWLSF(Long64_t, double, double, bool);
TH3D* threeDimMPPC(Long64_t, double, double, bool);
TH1D* totalLightYield();
TVector3 getMPPCPos(int);
void showAll(Long64_t, double, double);
bool passCosmicSelection(Long64_t, UInt_t, bool);



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
        data->GetMppc(randEvNum);
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
        data->GetMppc(randEvNum);
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
    vector<UInt_t> hitMPPC; 
    vector<UInt_t> hitADC; 

    while(true){
        randEvNum = rando->Rndm()*entries;
        if( passCosmicSelection(randEvNum, lyLim, true) ) break;
    }
    showAll(randEvNum, 1, 1000);
}


void showAll(Long64_t evNum, double limLow=1, double limHigh=300){

    vector<TH2D*> listOfHist = twoFaces(evNum, limLow, limHigh, true);
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



vector<TH2D*> twoFaces(Long64_t evNum, double limLow=1, double limHigh=300, bool noDraw=false){

    //Setup hhistos
    TH2D* facexzHisto = new TH2D("xz", "XZ View;X [cm]; Z [cm]", 7,0,7, 7,0,7);
    TH2D* faceyzHisto = new TH2D("yz", "YZ View;Y [cm]; Z [cm]", 7,0,7, 7,0,7);

    facexzHisto->GetZaxis()->SetRangeUser(limLow, limHigh);
    faceyzHisto->GetZaxis()->SetRangeUser(limLow, limHigh);

	Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    data->GetMppc(evNum);
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

    vector<TH2D*> listOfHist; 
    listOfHist.push_back(facexzHisto);
    listOfHist.push_back(faceyzHisto);

    return(listOfHist);

}

TH3D* threeDimWLSF(Long64_t evNum, double limLow=1, double limHigh=300, bool noDraw=false){
    //Setup hhistos
    TH3D* histo = new TH3D("xyzWLSF", "3D WLSF View;X [mm]; Y [mm]; Z[mm]", 14,0,70, 14,0,70 , 14,0,70);

    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    data->GetMppc(evNum);
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
    data->GetMppc(evNum);
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

    TH1D* lyHisto = new TH1D("totalLightYield", "totalLightYield", 100, 0, 6000);
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
        double totalLY = 0;
        for (int i=0;i<18;i++) totalLY += data->ADC(i);
        lyHisto->Fill(totalLY);
	}
    TCanvas* c1 = new TCanvas();
    lyHisto->Draw();

    return lyHisto;
}

vector<TH1D*> perMPPCLightYield()
{
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();

    vector<TH1D*> listOfHist; 
    for (int i=1;i<=18;i++){
        TH1D* lyHisto = new TH1D(Form("lyHisto%d",i), Form("MPPC %d; Light Yield [ADC]; Counts",i), 300, 0, 6000);
        listOfHist.push_back(lyHisto);
    }

    for( Long64_t ientry=0; ientry<entries; ++ientry )
    {
        data->GetMppc(ientry);
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

vector<TH1D*> perMPPCCosmicLightYieldAndHits(UInt_t lyLim=1500)
{
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();

    TH1D* hitsHisto = new TH1D("hitsHisto", "hitsHisto;MPPC Number; Hits From Cosmics", 18, 1, 19);

    vector<TH1D*> listOfHist; 
    for (int i=1;i<=18;i++){
        TH1D* lyHisto = new TH1D(Form("lyHisto%d",i), Form("MPPC %d; Light Yield [ADC]; Counts",i), 300, lyLim, 6000);
        listOfHist.push_back(lyHisto);
    }

    for( Long64_t ientry=0; ientry<entries; ++ientry )
    {
        // Check if we pass the cosmic selection for this event
        if( !passCosmicSelection(ientry, lyLim, false) ) continue;

        data->GetMppc(ientry);
        for (int i=0;i<18;i++){
            if(data->ADC(i) > lyLim){
                gPad->SetLogy();
                hitsHisto->Fill(i+1);
                listOfHist[i]->Fill(data->ADC(i));
            }
        }
    }

    //Setup histos
    TH2D* facexzHisto = new TH2D("MeanCosmicLY_XZ", "Mean Cosmic LY;X [cm]; Z [cm]", 7,0,7, 7,0,7);
    TH2D* faceyzHisto = new TH2D("MeanCosmicLY_YZ", "Mean Cosmic LY;Y [cm]; Z [cm]", 7,0,7, 7,0,7);

    facexzHisto->GetZaxis()->SetRangeUser(lyLim, 3000);
    faceyzHisto->GetZaxis()->SetRangeUser(lyLim, 3000);

    //Loop over the MPPCs 
    for(int i=1; i<=18; i++){
        TVector3 pos = getMPPCPos(i);
        //Face 1 or Face 3
        if(pos.Y()==0 || pos.Y()==7){
            facexzHisto->Fill(pos.X(), pos.Z(), listOfHist[i-1]->GetMean());
        }
        //Face 2 or Face 4
        else if(pos.X()==0 || pos.X()==7){
            faceyzHisto->Fill(pos.Y(), pos.Z(), listOfHist[i-1]->GetMean());
        }
        else{
            std::cout << "WARNING: found invalid MPPC position for MPPC " << i << std::endl;
            pos.Print();
        }
    }

    gStyle->SetOptTitle(1); 
    TCanvas *c1 = new TCanvas("c1","c1",1900,950);
    gStyle->SetOptStat(0);
    c1->Divide(6,3);
    for (int i=0;i<18;i++){
        c1->cd(i+1);
        listOfHist[i]->Draw();
    }
    TCanvas *c2 = new TCanvas();
    hitsHisto->Draw();

    TCanvas *c3 = new TCanvas("c3","multipads",1600,700);
    gStyle->SetOptStat(0);
    c3->Divide(2,1);
    c3->cd(1);
    facexzHisto->Draw("colzTEXT45");
    c3->cd(2);
    faceyzHisto->Draw("colzTEXT45");

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
        data->GetMppc(ientry);
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


bool passCosmicSelection(Long64_t evNum, UInt_t lyLim, bool verbose){
    bool pass = false;

    Mppc *data = filereader.tmCD(); 
    data->GetMppc(evNum);

    UInt_t totalLY = 0;
    UInt_t nHits = 0;
    UInt_t nHitsTopX = 0;
    UInt_t nHitsTopY = 0;
    UInt_t nHitsBotX = 0;
    UInt_t nHitsBotY = 0;
    bool isHit = false;
    vector<UInt_t> hitMPPC; 
    vector<UInt_t> hitADC; 

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
    //if(nHits>=6 && nHitsTopX==1 && nHitsTopY==1 && nHitsBotX==1 && nHitsBotY==1) pass = true;
    if(nHitsTopX==1 && nHitsTopY==1 && nHitsBotX==1 && nHitsBotY==1) pass = true;
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
