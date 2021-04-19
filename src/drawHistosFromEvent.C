R__LOAD_LIBRARY(fromUmut/GluedCubes/Analysis/TreeManager_C.so);

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

TreeManager filereader("./fromUmut/GluedCubes/SpecialRun_GluedCubesNoTeflon_15April2021.root");

using namespace std;

//gStyle->SetPaintTextFormat("6.1f");

//Function prototypes
vector<TH2D*> twoFaces(Long64_t, double, double, bool);
TH3D* threeDimWLSF(Long64_t, double, double, bool);
TH3D* threeDimMPPC(Long64_t, double, double, bool);
TH1D* totalLightYield();
TVector3 getMMPCPos(int);
void showAll(Long64_t, double, double);


void showHighLYEvent(double lyLim, int seed=0){
    //Initalise RNG generator
    TRandom3* rando = new TRandom3();
    rando->SetSeed(seed);
    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    int randEvNum = 0;

    double totalLY = 0;
    while(true){
        randEvNum = rando->Rndm()*entries;
        data->GetMppc(randEvNum);
        totalLY = 0;
        for (int i=0;i<18;i++) totalLY += data->ADC(i);
        if(totalLY>lyLim) break;
    }

    std::cout << "Chosen event number: " << randEvNum << std::endl;
    std::cout << "  Light Yield (ADC): " << totalLY << std::endl;

    showAll(randEvNum, 1, 500);
}


void showAll(Long64_t evNum, double limLow=1, double limHigh=300){

    vector<TH2D*> listOfHist = twoFaces(evNum, limLow, limHigh, true);
    TH3D* histoWLSF = threeDimWLSF(evNum, limLow, limHigh, true);
    TH3D* histoMPPC = threeDimMPPC(evNum, limLow, limHigh, true);

    TCanvas *c1 = new TCanvas("c1","multipads",950,950);
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
    TH2D* facexzHisto = new TH2D("xz", "xz;X [cm]; Z [cm]", 7,0,7, 7,0,7);
    TH2D* faceyzHisto = new TH2D("yz", "yz;Y [cm]; Z [cm]", 7,0,7, 7,0,7);

    facexzHisto->GetZaxis()->SetRangeUser(limLow, limHigh);
    faceyzHisto->GetZaxis()->SetRangeUser(limLow, limHigh);

	Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    data->GetMppc(evNum);
    std::cout << "Loading event number: " << evNum << std::endl;

    //Loop over the MPPCs in the event
    for(int i=1; i<=18; i++){
        TVector3 pos = getMMPCPos(i);
        //Face 1 or Face 3
        if(pos.Y()==0 || pos.Y()==7){
            facexzHisto->Fill(pos.X(), pos.Z(), data->ADC(i));
        }
        //Face 2 or Face 4
        else if(pos.X()==0 || pos.X()==7){
            faceyzHisto->Fill(pos.Y(), pos.Z(), data->ADC(i));
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
    TH3D* histo = new TH3D("xyzWLSF", "xyzWLSF;X [mm]; Y [mm]; Z[mm]", 14,0,70, 14,0,70 , 14,0,70);

    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    data->GetMppc(evNum);
    std::cout << "Loading event number: " << evNum << std::endl;

    //Loop over the MPPCs in the event
    for(int i=1; i<=18; i++){
        TVector3 pos = getMMPCPos(i);
        //Face 1 or Face 3
        if(pos.Y()==0 || pos.Y()==7){
            for(int j=2; j<=12; j++){
              histo->Fill(10*pos.X()+5, 5*j, 10*pos.Z()+5, data->ADC(i));
            }
        }
        //Face 2 or Face 4
        else if(pos.X()==0 || pos.X()==7){
            for(int j=2; j<=12; j++){
              histo->Fill(5*j, 10*pos.Y()+5, 10*pos.Z()+5, data->ADC(i));
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
    TH3D* histo = new TH3D("xyzMPPC", "xyzMPPC;X [cm]; Y [cm]; Z[cm]", 7,0,7, 7,0,7 , 7,0,7);

    Mppc *data = filereader.tmCD(); 
    Long64_t entries = data->GetInputTree()->GetEntries();
    data->GetMppc(evNum);
    std::cout << "Loading event number: " << evNum << std::endl;

    //Loop over the MPPCs in the event
    for(int i=1; i<=18; i++){
        TVector3 pos = getMMPCPos(i);
        //Face 1 or Face 3
        if(pos.Y()==0 || pos.Y()==7){
            histo->Fill(pos.X(), 0.0, pos.Z(), data->ADC(i));
        }
        //Face 2 or Face 4
        else if(pos.X()==0 || pos.X()==7){
            histo->Fill(0.0, pos.Y(), pos.Z(), data->ADC(i));
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
}

TVector3 getMMPCPos(int mppcNum){
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
