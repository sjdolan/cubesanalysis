#include <iostream>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include "TCanvas.h"
#include "TTree.h"
#include "TROOT.h"
#include "TPostScript.h"
#include "TPDF.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TPad.h"
#include <vector>
#include<algorithm>

#include "TBox.h"
#include "TEveManager.h"
#include "TEveEventManager.h"
#include "TGLUtil.h"
#include "TEveViewer.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TGeoMaterial.h"
#include "TGeoManager.h"
#include "TSystem.h"
#include "TEveGeoNode.h"
//#include "TEveGeoTopNode.h"
#include "TEveBoxSet.h"
#include "TEveRGBAPalette.h"
#include "TGLViewer.h"
#include "TEveBrowser.h"
#include "TGTab.h"
#include "TColor.h"
#include "TKey.h"


using namespace std;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R__LOAD_LIBRARY(TreeManager_C.so);
//R__LOAD_LIBRARY(TreeManager_C.dll);
#include "TreeManager.h"
//gSystem->Load("TreeManager_C.so");
//gSystem->Load("TreeManager_C.dll");
//TreeManager filereader("mppc_cosmicmuon_firsttrial_29_March_2018.root");
//TreeManager filereader("mppc_cosmicmuon_secondtrial_29_March_2018.root");
//TreeManager filereader("mppc_cosmicmuon_or32_30_March_2018.root");
//TreeManager filereader("mppc_cosmicmuon_CH14-15toCH16-17_TH250_29_March_2018.root");
//TreeManager filereader("mppc_cosmicmuon_thirdtrial_30_March_2018.root");
//TreeManager filereader("data//scint_2-4_12-15.root");
//TreeManager filereader("data//scint_2-4_12-15.root");
//TreeManager filereader("data//scint_6-11.root");
//TreeManager filereader("data//scint_4-10.root");
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TreeManager filereader(_file0->GetName());
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#define DEBUG 0
#define DEBUG1 0

TCanvas *c  =0;
TCanvas *c1 =0;
TCanvas *c2 =0;
TCanvas *c3 =0;
TCanvas *c4 =0;
TCanvas *myCant =0;
TCanvas *myCants =0;
TCanvas *myCan =0;
TCanvas *myCanx =0;
TCanvas *c2b =0;
TCanvas *c2v =0;
TCanvas *c2vt =0;
TCanvas *c2vs =0;

TGNumberEntry *fNumberEntry74;
TGNumberEntry *fNumberEntry75;
TGNumberEntry *fNumberEntry755;
TGNumberEntry *fNumberEntry886;
TGNumberEntry *fNumberEntry8866;
TGNumberEntry *fNumberEntry8869;
//TGLabel *fLabel;
TGStatusBar *fStatusBar739;
TGRadioButton * fChanProbe;
TGCheckButton * fChannel1;
TGCheckButton * fChannel2;
TGCheckButton *fUpdateHisto;
TGCheckButton *fUpdateVCXO;

TGTextButton *fGain;
TGTextButton *fGainc;

TGLabel *fLabel7;
TGLabel *fLabel79;
TGLabel *fLabel79s;
Int_t EventNumber;
Int_t CurrentEventNumber;
Int_t RunNumber;
Int_t NumberOfRuns;
Int_t ChannelNumber;
Bool_t ApplyGainFitting;
Bool_t SelectedChannels;

std::vector<int> ConvertToBinary2nd(int entry, int number, int num_digits);
int ConvertToBinary(int entry, int number, int num_digits) ;

#if 0
struct myclass {
	bool operator() (double i,double j) { return (i<j);}
} myobject;

struct SelEvent{
	int x,y;
	double amp;
};
int cmp(const void *a, const void *b){
	struct SelEvent *a1 = (struct SelEvent *)a;
	struct SelEvent *a2 = (struct SelEvent *)b;
	cout << (*a1).amp << " " << (*a2).amp << endl;
	if((*a1).amp>(*a2).amp)return -1;
	else if((*a1).amp<(*a2).amp)return 1;
	else return 0;
};
class EventOrderByWAmp
{
public:
	bool operator()( const SelEvent& a, const SelEvent b ) const
	{
		return a.amp < b.amp;
	};
};
#endif

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
template <typename T>
void DisplayContents(const T& Input)
{
	for (auto iElement = Input.cbegin(); iElement != Input.cend(); ++iElement)
		cout << *iElement << endl;
	
	
}
*/

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



void FEB_Efficiency()
{
	cout << "Getting FEB efficiency" << endl;
	Mppc *data = filereader.tmCD();
	mppc->Process("FEB_Efficiency.C+");
}

void FEB_Behaviour()
{
	cout << "Getting FEB behaviour" << endl;
	Mppc *data = filereader.tmCD();
	mppc->Process("FEB_Behaviour.C+");
}

void FEB_Behaviour_WithPlots()
{
	cout << "Getting FEB behaviour with plots" << endl;
	Mppc *data = filereader.tmCD();
	mppc->Process("FEB_Behaviour_withplot.C+");
}


vector<int> GetNumberOfConnectedFEBs()
{
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	//
	set<int> MacFEB;
	int foo;
	vector<int> numFEB;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
		foo = data->Mac();
		MacFEB.insert(foo);
		//printf("mac5=0x%02x %d \n",data->Mac(), MacFEB);
	}
	for (set<int>::iterator it = MacFEB.begin() ; it != MacFEB.end(); ++it)
	numFEB.push_back(*it);
//	std::cout << ' ' << *it;
//	std::cout << '\n';
	cout << " ... " << numFEB.size() << endl;
	return numFEB;
}

vector<vector<int> > GetIndexes(int eve)
{
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	//
	vector<int> Index_FEB;
	vector<vector<int>> indx;
	int foo;
	int j;
	vector<int> febs = GetNumberOfConnectedFEBs();
	for(int i =0; i<febs.size(); i++)
	{
		indx.push_back(vector<int>());
		cout << "check for " << febs[i] << endl;
		j=0;
		for( Long64_t ientry=0; ientry<entries; ++ientry )
		{
			data->GetMppc(ientry);
			foo = data->Mac();
			if(foo==febs[i])
			{
				j++;
				if(j==eve)
				{
					indx.back().push_back(ientry);
					cout << eve << " " << foo << " " << ientry << endl;
				}
			}
		}
	}
	
	for(int i =0; i<indx.size(); i++)
	for(int k =0; k<indx.at(i).size(); k++)
		cout << " indx " << i << " indx[" << i << "][" << k << "]" << indx.at(i).at(k) << endl;
	return indx;
}	


void ShowEventOnCRT_CAEN_Multi(int ientry)
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	
	
	vector<int> febs = GetNumberOfConnectedFEBs();
	vector<vector<int> > indx = GetIndexes(ientry);
	c->Clear();
	c->Modified();
	c->Divide(2,3);
	c->Update();

	TH2I *hpl2d[10];
	for(int i =0; i<indx.size(); i++)
	{
		hpl2d[i] = new TH2I(Form("histEv_%i",i),"Cosmic muon event",10,-1,9,10,-1,9);
	}
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	int num, index = 0;
	Double_t ratio = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;

//
	for(int ii =0; ii<indx.size(); ii++)
	for(int k =0; k<indx.at(ii).size(); k++)
	{
		int sel = indx.at(ii).at(k);
		data->GetMppc(sel);
		std::vector<int> v_InCoincidence;
		TBox *Box[10];
		//
		
		index = 0;
		for(int i = 0; i < 8; i++)
		{
			for(int j = 8; j < 16; j++)
			{
				index++;
				BarX = data->ADC(2*i) + data->ADC(2*i+1);
				BarY = data->ADC(2*j) + data->ADC(2*j+1);
				if(BarX==0||BarY==0) continue;
				//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
				num = abs(8-j);
				if(BarX==0||BarY==0) 
				{
					ratio = 0;
				} else{
					ratio = (Double_t) BarX/BarY;
				}
				if(ratio>1.)
				{
					ratio = (Double_t) BarY/BarX;
				}
				hpl2d[ii]->SetBinContent(7-i+2,num+2,ratio*(BarX+BarY));
			}
		}
		
		for (int i = 0; i < 32; i++ )
		{
			if(i==31) v_InCoincidence = ConvertToBinary2nd(ientry,data->Coincidence(),16);
		}
		for(int i=0;i<v_InCoincidence.size();i++)
		{
			if(DEBUG1) cout << "flag " << v_InCoincidence[i] << endl;
			if(v_InCoincidence[i]<9)
				Box[i] = new TBox(9-v_InCoincidence[i],0,8-v_InCoincidence[i],8);
			if(v_InCoincidence[i]>8)
				Box[i] = new TBox(0,v_InCoincidence[i]-9,8,v_InCoincidence[i]-8);
		}
		
		c->cd(ii+1);
		hpl2d[ii]->GetYaxis()->SetBinLabel(9,"Scin16");
		hpl2d[ii]->GetYaxis()->SetBinLabel(8,"Scin15");
		hpl2d[ii]->GetYaxis()->SetBinLabel(7,"Scin14");
		hpl2d[ii]->GetYaxis()->SetBinLabel(6,"Scin13");
		hpl2d[ii]->GetYaxis()->SetBinLabel(5,"Scin12");
		hpl2d[ii]->GetYaxis()->SetBinLabel(4,"Scin11");
		hpl2d[ii]->GetYaxis()->SetBinLabel(3,"Scin10");
		hpl2d[ii]->GetYaxis()->SetBinLabel(2,"Scin9");
		hpl2d[ii]->GetYaxis()->SetLabelSize(0.04);
		//
		hpl2d[ii]->GetXaxis()->SetBinLabel(9,"Scin1");
		hpl2d[ii]->GetXaxis()->SetBinLabel(8,"Scin2");
		hpl2d[ii]->GetXaxis()->SetBinLabel(7,"Scin3");
		hpl2d[ii]->GetXaxis()->SetBinLabel(6,"Scin4");
		hpl2d[ii]->GetXaxis()->SetBinLabel(5,"Scin5");
		hpl2d[ii]->GetXaxis()->SetBinLabel(4,"Scin6");
		hpl2d[ii]->GetXaxis()->SetBinLabel(3,"Scin7");
		hpl2d[ii]->GetXaxis()->SetBinLabel(2,"Scin8");
		hpl2d[ii]->GetXaxis()->SetLabelSize(0.04);
		hpl2d[ii]->SetStats(kFALSE);
		
		hpl2d[ii]->SetTitle(Form("Cosmic Event %i in FEB_Mac: %i",ientry,febs[ii]));
		hpl2d[ii]->GetXaxis()->SetTitleOffset(1.5);
		hpl2d[ii]->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
		hpl2d[ii]->GetYaxis()->SetTitleOffset(1.5);
		hpl2d[ii]->GetYaxis()->SetTitle("Top Layer (10 mm)");
		hpl2d[ii]->DrawCopy("col text");
		for(int i=0;i<v_InCoincidence.size();i++)
		{
			Box[i]->SetFillStyle(0);
			Box[i]->SetLineColor(2);
			Box[i]->SetLineWidth(2);
			Box[i]->Draw("same");
		}
		//gPad->Modified();
		//gPad->Update();
		c->Update();
		c->Modified();
		c->Update();
		hpl2d[ii]->GetZaxis()->SetRangeUser(min, max);
		
	}
//	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




void GetRun2_4vs12_15()
{
	cout << "Getting scatter plot for Run2_4vs12_15" << endl;
	TCanvas *myCan = new TCanvas("Canvas","Quality Check",200,10,700,500);
	Mppc *data = filereader.tmCD();
	myCan->Divide(2,2);
	myCan->Update();
	myCan->Modified();
	TH2F *hist2[4];
	myCan->cd(1);
	mppc->Draw("chg[2]+chg[3]:chg[22]+chg[23]>> htemp(100,0,8200,100,0,8200)");
	hist2[0] = (TH2F*)gDirectory->Get("htemp");
	hist2[0]->SetTitle("Scin #2 [ch2+ch3] vs Scin #12 [ch22+ch23];Scin 12 [ch22+ch23];Scin 2 [ch2+ch3]");
	hist2[0]->SetMarkerStyle(20);
	hist2[0]->SetMarkerSize(0.2);
	hist2[0]->SetStats(kFALSE);
	hist2[0]->Draw();
	hist2[0]->GetXaxis()->SetTitleOffset(1.4); 
	hist2[0]->GetYaxis()->SetTitleOffset(1.4); 
	myCan->Update();
	myCan->cd(2);
	mppc->Draw("chg[2]+chg[3]:chg[28]+chg[29]>> htemp2(100,0,8200,100,0,8200)");
	hist2[1] = (TH2F*)gDirectory->Get("htemp2");
	hist2[1]->SetTitle("Scin #2 [ch2+ch3] vs Scin #15 [ch22+ch23];Scin 15 [ch28+ch29];Scin 2 [ch2+ch3]");
	hist2[1]->SetMarkerStyle(20);
	hist2[1]->SetMarkerSize(0.2);
	hist2[1]->SetStats(kFALSE);
	hist2[1]->Draw();
	hist2[1]->GetXaxis()->SetTitleOffset(1.4); 
	hist2[1]->GetYaxis()->SetTitleOffset(1.4); 
	myCan->Update();
	myCan->cd(3);
	mppc->Draw("chg[6]+chg[7]:chg[22]+chg[23]>> htemp3(100,0,8200,100,0,8200)");
	hist2[2] = (TH2F*)gDirectory->Get("htemp3");
	hist2[2]->SetTitle("Scin #4 [ch6+ch7] vs Scin #12 [ch22+ch23];Scin 12 [ch23+ch24];Scin 4 [ch6+ch7]");
	hist2[2]->SetMarkerStyle(20);
	hist2[2]->SetMarkerSize(0.2);
	hist2[2]->SetStats(kFALSE);
	hist2[2]->Draw();
	hist2[2]->GetXaxis()->SetTitleOffset(1.4); 
	hist2[2]->GetYaxis()->SetTitleOffset(1.4); 
	myCan->Update();
	myCan->cd(4);
	mppc->Draw("chg[6]+chg[7]:chg[28]+chg[29]>> htemp4(100,0,8200,100,0,8200)");
	hist2[3] = (TH2F*)gDirectory->Get("htemp4");
	hist2[3]->SetTitle("Scin #4 [ch6+ch7] vs Scin #15 [ch28+ch29];Scin 15 [ch28+ch29];Scin 4 [ch6+ch7]");
	hist2[3]->SetMarkerStyle(20);
	hist2[3]->SetMarkerSize(0.2);
	hist2[3]->SetStats(kFALSE);
	hist2[3]->Draw();
	hist2[3]->GetXaxis()->SetTitleOffset(1.4); 
	hist2[3]->GetYaxis()->SetTitleOffset(1.4); 
	myCan->Update();
	myCan->SaveAs("Scin2_4vs12_15.ps");
}
/////////////////////////////////////////////

TH2F *hist2d[16];
void TwoFibersInBar1()
{

  	for(int i = 0; i<16; i++)
	{
		TH2F *m = (TH2F*)gROOT->FindObject(Form("hist2d_%i",i));
		if(m) m->Delete();
	}
	
	for(int i=0;i<16;i++)
	{
		const char* name = Form("hist2d_%i",i);
		hist2d[i] = new TH2F(name,name,1000,0,4050,1000,0,4050);
	}
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	cout << "Number of entries " << endl;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
		//
		for(int i = 0; i < 16; i++)
		{
			hist2d[i]->Fill(data->ADC(2*i),data->ADC(2*i+1));
		}
	}
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Quality Check",200,10,700,500);
	for(int i=0; i < 16; i++)
	{
		hist2d[i]->SetLineWidth(2);
		hist2d[i]->SetMarkerStyle(2);
		hist2d[i]->SetStats(kFALSE);
		hist2d[i]->SetTitle(Form("Two fiber in Bar #%i;Right Fiber [CH%i];Left Fiber [CH%i]",i+1,2*i,2*i+1));
		hist2d[i]->Draw();
		if(i == 0 ) myCan->SaveAs("TwoFibersInBar.pdf(");
		if(i == 15 ) myCan->SaveAs("TwoFibersInBar.pdf)");
		if(i > 0 && i < 15) myCan->SaveAs("TwoFibersInBar.pdf");

	}
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void TwoLayersInModule()
{

  	for(int i = 0; i<8; i++)
	{
		TH2F *m = (TH2F*)gROOT->FindObject(Form("hpl2d_%i",i));
		if(m) m->Delete();
	}

	TH2F *hpl2d[8][8];
	for(int i = 0; i < 8; i++)
	{
		for(int j = 0; j < 8; j++)
		{
			const char* name = Form("hpl2d_%i_%i",i,j);
			hpl2d[i][j] = new TH2F(name,name,1000,0,8100,1000,0,8100);
		}
	}
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	cout << "Number of entries " << entries << endl;
	UInt_t BarX, BarY;
	int num, index = 0;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
		//
		for(int i = 0; i < 8; i++)
		{
			for(int j = 8; j < 16; j++)
			{
				index++;
				BarX = data->ADC(2*i) + data->ADC(2*i+1);
				BarY = data->ADC(2*j) + data->ADC(2*j+1);
				//cout << 2*i << " " << 2*i+1 << " " << 2
				num = abs(8-j);
				hpl2d[i][num]->Fill(BarX,BarY);
				if(BarX>200&&BarY>200)
				cout << ientry << " " << i << " " << num << " " << BarX << " " << BarY << endl;
			}
		}
	}

	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Quality Check",200,10,700,500);
	int nn=0;
	for(int i = 0; i < 8; i++)
	{
		for(int j = 0; j < 8; j++)
		{
			nn++;
			hpl2d[i][j]->SetLineWidth(2);
			hpl2d[i][j]->SetMarkerStyle(2);
			hpl2d[i][j]->SetStats(kFALSE);
			hpl2d[i][j]->SetTitle(Form("15 mm Bar-%i vs 10 mm Bar-%i;15 mm Bar [CH%i+CH%i];10 mm Bar [CH%i+CH%i]",i+1,j+1,2*i,2*i+1,2*j,2*j+1));
			hpl2d[i][j]->Draw();
			if(nn==1) myCan->SaveAs("15mmBarTo10mmBar.pdf(");
			if(nn == index) myCan->SaveAs("15mmBarTo10mmBar.pdf)");
			if(nn>1&&nn<index) myCan->SaveAs("15mmBarTo10mmBar.pdf");
		}
	}
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ShowEventOnCRT(int ientry){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
//	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
//	if (C) delete C;
	
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	
	//TCanvas *c = new TCanvas("Canvas","Show Event on CRT",200,10,700,500);
	//TPostScript *ps = new TPostScript(Form("CRT_Event_%i.ps",ientry),112);

	TH2I *hpl2d = new TH2I("histEv","Cosmic muon event",10,-1,9,10,-1,9);

	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	int num, index = 0;
	Double_t ratio = 0.;
	Double_t ratio_t = 0.;
	Double_t ratio_b = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;
	data->GetMppc(ientry);
	//
	for(int i = 0; i < 8; i++)
	{
		for(int j = 8; j < 16; j++)
		{
			index++;
			BarX = data->ADC(2*i) + data->ADC(2*i+1);
			BarY = data->ADC(2*j) + data->ADC(2*j+1);
			//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
			num = abs(8-j);
			if(BarX==0||BarY==0) 
			{
				ratio = 0;
			} else{
				ratio = (Double_t) BarX/BarY;
			}
			if(ratio>1.)
			{
				ratio = (Double_t) BarY/BarX;
			}
			//hpl2d->SetBinContent(i+2,num+2,ratio*(BarX+BarY));
			hpl2d->SetBinContent(7-i+2,num+2,ratio*(BarX+BarY));
		}
	}
	//ps->NewPage();
	c->cd();
	hpl2d->GetYaxis()->SetBinLabel(9,"Scin16");
	hpl2d->GetYaxis()->SetBinLabel(8,"Scin15");
	hpl2d->GetYaxis()->SetBinLabel(7,"Scin14");
	hpl2d->GetYaxis()->SetBinLabel(6,"Scin13");
	hpl2d->GetYaxis()->SetBinLabel(5,"Scin12");
	hpl2d->GetYaxis()->SetBinLabel(4,"Scin11");
	hpl2d->GetYaxis()->SetBinLabel(3,"Scin10");
	hpl2d->GetYaxis()->SetBinLabel(2,"Scin9");
	hpl2d->GetYaxis()->SetLabelSize(0.04);
	//
	hpl2d->GetXaxis()->SetBinLabel(9,"Scin1");
	hpl2d->GetXaxis()->SetBinLabel(8,"Scin2");
	hpl2d->GetXaxis()->SetBinLabel(7,"Scin3");
	hpl2d->GetXaxis()->SetBinLabel(6,"Scin4");
	hpl2d->GetXaxis()->SetBinLabel(5,"Scin5");
	hpl2d->GetXaxis()->SetBinLabel(4,"Scin6");
	hpl2d->GetXaxis()->SetBinLabel(3,"Scin7");
	hpl2d->GetXaxis()->SetBinLabel(2,"Scin8");
	hpl2d->GetXaxis()->SetLabelSize(0.04);
	hpl2d->SetStats(kFALSE);

	hpl2d->SetTitle(Form("Cosmic Event %i",ientry));
	hpl2d->GetXaxis()->SetTitleOffset(1.5);
	hpl2d->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
	hpl2d->GetYaxis()->SetTitleOffset(1.5);
	hpl2d->GetYaxis()->SetTitle("Top Layer (10 mm)");
	hpl2d->DrawCopy("col text");
	//gPad->Modified();
	//gPad->Update();
	hpl2d->GetZaxis()->SetRangeUser(min, max);
	c->Update();
	c->Modified();
	c->Update();
	//
	
	//ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ShowEventOnCRT_cut(int ientry){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	

	TCanvas *myCan = new TCanvas("Canvas","Show Event on CRT",200,10,700,500);
	TPostScript *ps = new TPostScript(Form("CRT_Event_%i.ps",ientry),112);

	TH2I *hpl2d = new TH2I("histEv","Cosmic muon event",10,-1,9,10,-1,9);
	Int_t  nbinsx = hpl2d->GetXaxis()->GetNbins();
	Int_t  nbinsy = hpl2d->GetYaxis()->GetNbins();
	Int_t ibin,bin;

	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	UInt_t adc_value;
	int count1 = 0, count2 = 0, count3 = 0;
	int num, index = 0;
	Double_t ratio = 0.;
	Double_t ratio_t = 0.;
	Double_t ratio_b = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;
	data->GetMppc(ientry);
	//
	for(int i = 0; i < 8; i++)
	{
		for(int j = 8; j < 16; j++)
		{
			index++;
			BarX = data->ADC(2*i) + data->ADC(2*i+1);
			BarY = data->ADC(2*j) + data->ADC(2*j+1);
			//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
			num = abs(8-j);
			if(BarX==0||BarY==0) 
			{
				ratio = 0;
			} else{
				ratio = (Double_t) BarX/BarY;
			}
			if(ratio>1.)
			{
				ratio = (Double_t) BarY/BarX;
			}
			//hpl2d->SetBinContent(i+2,num+2,ratio*(BarX+BarY));
			hpl2d->SetBinContent(7-i+2,num+2,ratio*(BarX+BarY));
		}
	}
	ps->NewPage();
	hpl2d->GetYaxis()->SetBinLabel(9,"Scin16");
	hpl2d->GetYaxis()->SetBinLabel(8,"Scin15");
	hpl2d->GetYaxis()->SetBinLabel(7,"Scin14");
	hpl2d->GetYaxis()->SetBinLabel(6,"Scin13");
	hpl2d->GetYaxis()->SetBinLabel(5,"Scin12");
	hpl2d->GetYaxis()->SetBinLabel(4,"Scin11");
	hpl2d->GetYaxis()->SetBinLabel(3,"Scin10");
	hpl2d->GetYaxis()->SetBinLabel(2,"Scin9");
	hpl2d->GetYaxis()->SetLabelSize(0.04);
	//
	hpl2d->GetXaxis()->SetBinLabel(9,"Scin1");
	hpl2d->GetXaxis()->SetBinLabel(8,"Scin2");
	hpl2d->GetXaxis()->SetBinLabel(7,"Scin3");
	hpl2d->GetXaxis()->SetBinLabel(6,"Scin4");
	hpl2d->GetXaxis()->SetBinLabel(5,"Scin5");
	hpl2d->GetXaxis()->SetBinLabel(4,"Scin6");
	hpl2d->GetXaxis()->SetBinLabel(3,"Scin7");
	hpl2d->GetXaxis()->SetBinLabel(2,"Scin8");
	hpl2d->GetXaxis()->SetLabelSize(0.04);
	hpl2d->SetStats(kFALSE);
	
	hpl2d->SetTitle(Form("Cosmic Event %i",ientry));
	hpl2d->GetXaxis()->SetTitleOffset(1.5);
	hpl2d->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
	hpl2d->GetYaxis()->SetTitleOffset(1.5);
	hpl2d->GetYaxis()->SetTitle("Top Layer (10 mm)");
	hpl2d->DrawCopy("col text");
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hpl2d->GetZaxis()->SetRangeUser(min, max);
	ps->Close();
	
	
	for (Int_t biny=1;biny<=nbinsy;biny++)
	{
		for (Int_t binx=1;binx<=nbinsx;binx++)
		{
			ibin = hpl2d->GetBin(binx,biny);
			adc_value = hpl2d->GetBinContent(ibin);
			cout << "Bin[" << binx << "-" << biny << "]"
				 << " " << adc_value 
				 << endl;
			if(adc_value>1000) count1++;
			if(adc_value>1500) count2++;
			if(adc_value>2000) count3++;
		}
	}
	cout << "Number of hits with ADC_value " << endl;
	cout << "   >1000 --> " << count1 << endl;
	cout << "   >1500 --> " << count2 << endl;
	cout << "   >2000 --> " << count3 << endl;
	if(count1 > 1) cout << "Shower!!! " << ientry << endl;
	
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// following just for test
void ShowEventOnCRTTest(int ientry){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	

	TCanvas *myCan = new TCanvas("Canvas","Show Event on CRT",200,10,700,500);
	TPostScript *ps = new TPostScript(Form("CRT_Event_%i.ps",ientry),112);

	TH2I *hpl2d = new TH2I("histEv","Cosmic muon event",10,-1,9,10,-1,9);

	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	int num, index = 0;
	Double_t ratio = 0.;
	Double_t ratio_t = 0.;
	Double_t ratio_b = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;
	data->GetMppc(ientry);
	//
	for(int i = 0; i < 8; i++)
	{
		for(int j = 8; j < 16; j++)
		{
			index++;
			
			if(data->ADC(2*i)==0||data->ADC(2*i+1)==0)
			{ 
				ratio_b =0;
			}else{
				ratio_b = data->ADC(2*i) / data->ADC(2*i+1);
				if(ratio_b>1.)
					ratio_b = data->ADC(2*i+1) / data->ADC(2*i);
			}
			if(data->ADC(2*j)==0||data->ADC(2*j+1)==0)
			{ 
				ratio_t =0;
			}else{
				ratio_t = data->ADC(2*j) / data->ADC(2*j+1);
				if(ratio_t>1.)
					ratio_t = data->ADC(2*j+1) / data->ADC(2*j);
			}
			BarX = ratio_b*(data->ADC(2*i) + data->ADC(2*i+1));
			BarY = ratio_t*(data->ADC(2*j) + data->ADC(2*j+1));
			//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
			num = abs(8-j);
			if(BarX==0||BarY==0) 
			{
				ratio = 0;
			} else{
				ratio = (Double_t) BarX/BarY;
			}
			if(ratio>1.)
			{
				ratio = (Double_t) BarY/BarX;
			}
			//hpl2d->SetBinContent(i+2,num+2,ratio*(BarX+BarY));
			hpl2d->SetBinContent(7-i+2,num+2,ratio*(BarX+BarY));
		}
	}
	ps->NewPage();
	hpl2d->GetYaxis()->SetBinLabel(9,"Scin16");
	hpl2d->GetYaxis()->SetBinLabel(8,"Scin15");
	hpl2d->GetYaxis()->SetBinLabel(7,"Scin14");
	hpl2d->GetYaxis()->SetBinLabel(6,"Scin13");
	hpl2d->GetYaxis()->SetBinLabel(5,"Scin12");
	hpl2d->GetYaxis()->SetBinLabel(4,"Scin11");
	hpl2d->GetYaxis()->SetBinLabel(3,"Scin10");
	hpl2d->GetYaxis()->SetBinLabel(2,"Scin9");
	hpl2d->GetYaxis()->SetLabelSize(0.04);
	//
	hpl2d->GetXaxis()->SetBinLabel(9,"Scin1");
	hpl2d->GetXaxis()->SetBinLabel(8,"Scin2");
	hpl2d->GetXaxis()->SetBinLabel(7,"Scin3");
	hpl2d->GetXaxis()->SetBinLabel(6,"Scin4");
	hpl2d->GetXaxis()->SetBinLabel(5,"Scin5");
	hpl2d->GetXaxis()->SetBinLabel(4,"Scin6");
	hpl2d->GetXaxis()->SetBinLabel(3,"Scin7");
	hpl2d->GetXaxis()->SetBinLabel(2,"Scin8");
	hpl2d->GetXaxis()->SetLabelSize(0.04);
	hpl2d->SetStats(kFALSE);
	
	hpl2d->SetTitle(Form("Cosmic Event %i",ientry));
	hpl2d->GetXaxis()->SetTitleOffset(1.5);
	hpl2d->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
	hpl2d->GetYaxis()->SetTitleOffset(1.5);
	hpl2d->GetYaxis()->SetTitle("Top Layer (10 mm)");
	hpl2d->DrawCopy("col text");
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hpl2d->GetZaxis()->SetRangeUser(min, max);
	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ShowEventOnCRT_CAEN(int ientry)
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
//	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
//	if (C) delete C;
	
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	
//	TCanvas *myCan = new TCanvas("Canvas","Show Event on CRT",200,10,700,500);
//	TPostScript *ps = new TPostScript(Form("CRT_Event_%i.ps",ientry),112);

	TH2I *hpl2d = new TH2I("histEv","Cosmic muon event",10,-1,9,10,-1,9);

	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	int num, index = 0;
	Double_t ratio = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;
	data->GetMppc(ientry);
	std::vector<int> v_InCoincidence;
	TBox *Box[10];
	//
	for(int i = 0; i < 8; i++)
	{
		for(int j = 8; j < 16; j++)
		{
			index++;
			BarX = data->ADC(2*i) + data->ADC(2*i+1);
			BarY = data->ADC(2*j) + data->ADC(2*j+1);
			if(BarX==0||BarY==0) continue;
			//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
			num = abs(8-j);
			if(BarX==0||BarY==0) 
			{
				ratio = 0;
			} else{
				ratio = (Double_t) BarX/BarY;
			}
			if(ratio>1.)
			{
				ratio = (Double_t) BarY/BarX;
			}
			//hpl2d->SetBinContent(i+2,num+2,ratio*(BarX+BarY));
			hpl2d->SetBinContent(7-i+2,num+2,ratio*(BarX+BarY));
		}
	}

	for (int i = 0; i < 32; i++ )
	{
		if(i==31) v_InCoincidence = ConvertToBinary2nd(ientry,data->Coincidence(),16);
	}
//	cout << "number of flags: " << v_InCoincidence.size() << endl;
	for(int i=0;i<v_InCoincidence.size();i++)
	{
		if(DEBUG1) cout << "flag " << v_InCoincidence[i] << endl;
		if(v_InCoincidence[i]<9)
			Box[i] = new TBox(9-v_InCoincidence[i],0,8-v_InCoincidence[i],8);
		if(v_InCoincidence[i]>8)
			Box[i] = new TBox(0,v_InCoincidence[i]-9,8,v_InCoincidence[i]-8);
	}
	
	//ps->NewPage();
	c->cd();
	hpl2d->GetYaxis()->SetBinLabel(9,"Scin16");
	hpl2d->GetYaxis()->SetBinLabel(8,"Scin15");
	hpl2d->GetYaxis()->SetBinLabel(7,"Scin14");
	hpl2d->GetYaxis()->SetBinLabel(6,"Scin13");
	hpl2d->GetYaxis()->SetBinLabel(5,"Scin12");
	hpl2d->GetYaxis()->SetBinLabel(4,"Scin11");
	hpl2d->GetYaxis()->SetBinLabel(3,"Scin10");
	hpl2d->GetYaxis()->SetBinLabel(2,"Scin9");
	hpl2d->GetYaxis()->SetLabelSize(0.04);
	//
	hpl2d->GetXaxis()->SetBinLabel(9,"Scin1");
	hpl2d->GetXaxis()->SetBinLabel(8,"Scin2");
	hpl2d->GetXaxis()->SetBinLabel(7,"Scin3");
	hpl2d->GetXaxis()->SetBinLabel(6,"Scin4");
	hpl2d->GetXaxis()->SetBinLabel(5,"Scin5");
	hpl2d->GetXaxis()->SetBinLabel(4,"Scin6");
	hpl2d->GetXaxis()->SetBinLabel(3,"Scin7");
	hpl2d->GetXaxis()->SetBinLabel(2,"Scin8");
	hpl2d->GetXaxis()->SetLabelSize(0.04);
	hpl2d->SetStats(kFALSE);
	
	hpl2d->SetTitle(Form("Cosmic Event %i",ientry));
	hpl2d->GetXaxis()->SetTitleOffset(1.5);
	hpl2d->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
	hpl2d->GetYaxis()->SetTitleOffset(1.5);
	hpl2d->GetYaxis()->SetTitle("Top Layer (10 mm)");
	hpl2d->DrawCopy("col text");
	for(int i=0;i<v_InCoincidence.size();i++)
	{
		Box[i]->SetFillStyle(0);
		Box[i]->SetLineColor(2);
		Box[i]->SetLineWidth(2);
		Box[i]->Draw("same");
	}
	//gPad->Modified();
	//gPad->Update();
	c->Update();
	c->Modified();
	c->Update();
	hpl2d->GetZaxis()->SetRangeUser(min, max);
//	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void CRTMuonHits()
{
	cout << "Cosmic Muon Hits on CRT" << endl;
	Mppc *data = filereader.tmCD();
	mppc->Process("FEB_CRTScatter.C+");
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ShowAllEventOnCRT()
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	

	TCanvas *myCan = new TCanvas("Canvas","Show Event on CRT",200,10,700,500);
	TPostScript *ps = new TPostScript("CRT_AllEvent.ps",112);

	TH2I *hpl2d = new TH2I("histEv","Cosmic muon event",10,-1,9,10,-1,9);

	Int_t  nbinsx = hpl2d->GetXaxis()->GetNbins();
	Int_t  nbinsy = hpl2d->GetYaxis()->GetNbins();
	Int_t ibin,bin;
	for (Int_t biny=1;biny<=nbinsy;biny++) {
		for (Int_t binx=1;binx<=nbinsx;binx++) {
			ibin = hpl2d->GetBin(binx,biny);
			hpl2d->SetBinContent(ibin,0);
		}
	}
	//
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	int num;
	Double_t ratio = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;
	std::vector<int> v_InCoincidence;
	std::vector<double> sel;
	TBox *Box[10];
	//
	Double_t Sel[64];
	Int_t SelIndex[64];
	Int_t DirectionX[64];
	Int_t DirectionY[64];
	Int_t size;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
		//if(ientry>100) break;
		size = 0;
		for(int i = 0; i < 8; i++)
		{
			for(int j = 8; j < 16; j++)
			{
				BarX = data->ADC(2*i) + data->ADC(2*i+1);
				BarY = data->ADC(2*j) + data->ADC(2*j+1);
				if(BarX==0||BarY==0) continue;
				//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
				num = abs(8-j);
				ratio = (Double_t) BarX/BarY;
				if(ratio>1.)
				{
					ratio = (Double_t) BarY/BarX;
				}
				if( DEBUG ) 
					cout << "---> " << ratio*(BarX+BarY) << " " << i << " " << num << " " << size << " :: " << 8*i+num<< endl;
				Sel[size] = ratio*(BarX+BarY);
				DirectionX[size] = i;
				DirectionY[size] = num;
				SelIndex[size] = size;
				size++;
			}
		}
		bool switched = false;
		do
		{
			switched = false;
			for(int i = 1; i < size; i++)
			{
				if(Sel[SelIndex[i - 1]] < Sel[SelIndex[i]])
				{
					int temp = SelIndex[i];
					SelIndex[i] = SelIndex[i - 1];
					SelIndex[i - 1] = temp;
					switched = true;
				}
			}
		}
		while(switched);
		for(int i=0;i<size;i++)
		{
			if( DEBUG )
				cout << "sorted " << SelIndex[i] << " " << Sel[SelIndex[i]] << "--> " << DirectionX[SelIndex[i]] << " " << DirectionY[SelIndex[i]] << endl;
		}
		bin = hpl2d->GetBin(7-DirectionX[SelIndex[0]]+2,DirectionY[SelIndex[0]]+2);
		// The maximum value is the first sorted element
		//hpl2d->SetBinContent(7-DirectionX[SelIndex[0]]+2,DirectionY[SelIndex[0]]+2,Sel[SelIndex[0]]);
		//hpl2d->SetBinContent(7-DirectionX[SelIndex[0]]+2,DirectionY[SelIndex[0]]+2,1);
		hpl2d->AddBinContent(bin,1);

	}
	ps->NewPage();
	hpl2d->GetYaxis()->SetBinLabel(9,"Scin16");
	hpl2d->GetYaxis()->SetBinLabel(8,"Scin15");
	hpl2d->GetYaxis()->SetBinLabel(7,"Scin14");
	hpl2d->GetYaxis()->SetBinLabel(6,"Scin13");
	hpl2d->GetYaxis()->SetBinLabel(5,"Scin12");
	hpl2d->GetYaxis()->SetBinLabel(4,"Scin11");
	hpl2d->GetYaxis()->SetBinLabel(3,"Scin10");
	hpl2d->GetYaxis()->SetBinLabel(2,"Scin9");
	hpl2d->GetYaxis()->SetLabelSize(0.04);
	//
	hpl2d->GetXaxis()->SetBinLabel(9,"Scin1");
	hpl2d->GetXaxis()->SetBinLabel(8,"Scin2");
	hpl2d->GetXaxis()->SetBinLabel(7,"Scin3");
	hpl2d->GetXaxis()->SetBinLabel(6,"Scin4");
	hpl2d->GetXaxis()->SetBinLabel(5,"Scin5");
	hpl2d->GetXaxis()->SetBinLabel(4,"Scin6");
	hpl2d->GetXaxis()->SetBinLabel(3,"Scin7");
	hpl2d->GetXaxis()->SetBinLabel(2,"Scin8");
	hpl2d->GetXaxis()->SetLabelSize(0.04);
	hpl2d->SetStats(kFALSE);
	
	//hpl2d->SetTitle(Form("Cosmic Event %i",ientry));
	hpl2d->GetXaxis()->SetTitleOffset(1.5);
	hpl2d->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
	hpl2d->GetYaxis()->SetTitleOffset(1.5);
	hpl2d->GetYaxis()->SetTitle("Top Layer (10 mm)");
	hpl2d->DrawCopy("col text");
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hpl2d->GetZaxis()->SetRangeUser(min, max);
	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ShowAllEventOnCRT_Fraction()
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	TH2I *mt = (TH2I*) gROOT->FindObject("histEvTop");
	if(mt) delete mt;
	TH2I *mb = (TH2I*) gROOT->FindObject("histEvBot");
	if(mb) delete mb;
	

	TCanvas *myCan = new TCanvas("Canvas","Show Event on CRT",200,10,700,500);
	TPostScript *ps = new TPostScript("CRT_AllEvent_Fraction.ps",112);

	TH2I *hpl2d = new TH2I("histEv","Cosmic muon event",10,-1,9,10,-1,9);
	TH2I *hpltop = new TH2I("histTop","Cosmic muon event (Top Layer Signal Arrives First)",10,-1,9,10,-1,9);
	TH2I *hplbot = new TH2I("histBot","Cosmic muon event (Bottom Layer Signal Arrives First)",10,-1,9,10,-1,9);
	TH2I *hpleq = new TH2I("histEq","Cosmic muon event (Both Signal has same strength)",10,-1,9,10,-1,9);

	Int_t  nbinsx = hpl2d->GetXaxis()->GetNbins();
	Int_t  nbinsy = hpl2d->GetYaxis()->GetNbins();
	Int_t ibin,bin,ibint,bint,ibinb,binb,ibineq,bineq;
	for (Int_t biny=1;biny<=nbinsy;biny++) {
		for (Int_t binx=1;binx<=nbinsx;binx++) {
			ibin = hpl2d->GetBin(binx,biny);
			hpl2d->SetBinContent(ibin,0);
			ibint = hpltop->GetBin(binx,biny);
			hpltop->SetBinContent(ibint,0);
			ibinb = hplbot->GetBin(binx,biny);
			hplbot->SetBinContent(ibinb,0);
			ibineq = hpleq->GetBin(binx,biny);
			hpleq->SetBinContent(ibineq,0);
		}
	}
	//
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	int num;
	Double_t ratio = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;
	std::vector<int> v_InCoincidence;
	std::vector<double> sel;
	TBox *Box[10];
	UInt_t PossibleCoincidence;
	UInt_t PossibleCoincidence2;
	Double_t ratioBottom[8];
	Double_t weightedBottom[8];
	Double_t ratioTop[8];
	Double_t weightedTop[8];
	//
	Double_t Sel[64];
	Int_t SelIndex[64];
	Int_t DirectionX[64];
	Int_t DirectionY[64];
	Int_t size;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
		if(ientry>1000) break;
		size = 0;
		for(int i = 0; i < 8; i++)
		{
			if (data->ADC(2*i) == 0 || data->ADC(2*i+1) ==0)
			{
				ratioBottom[i]    = 0;
				weightedBottom[i] = 0;
			}else{
				ratioBottom[i] = (Double_t)data->ADC(2*i)/data->ADC(2*i+1);
				if(ratioBottom[i]>1.) 
				{
					ratioBottom[i] = (Double_t)data->ADC(2*i+1)/data->ADC(2*i);
				}
				weightedBottom[i] = ratioBottom[i]*(data->ADC(2*i)+data->ADC(2*i+1));
			}
			if (data->ADC(2*i+16) == 0 || data->ADC(2*i+1+16) ==0)
			{
				ratioTop[i]    = 0;
				weightedTop[i] = 0;
			}else{
				ratioTop[i] = (Double_t)data->ADC(2*i+16)/data->ADC(2*i+1+16);
				if(ratioTop[i]>1.) 
				{
					ratioTop[i] = (Double_t)data->ADC(2*i+1+16)/data->ADC(2*i+16);
				}
				weightedTop[i] = ratioTop[i]*(data->ADC(2*i+16)+data->ADC(2*i+1+16));
			}
			///////
			for(int j = 8; j < 16; j++)
			{
				BarX = data->ADC(2*i) + data->ADC(2*i+1);
				BarY = data->ADC(2*j) + data->ADC(2*j+1);
				if(BarX==0||BarY==0) continue;
				//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
				num = abs(8-j);
				ratio = (Double_t) BarX/BarY;
				if(ratio>1.)
				{
					ratio = (Double_t) BarY/BarX;
				}
				if( DEBUG ) 
					cout << "---> " << ratio*(BarX+BarY) << " " << i << " " << num << " " << size << " :: " << 8*i+num<< endl;
				Sel[size] = ratio*(BarX+BarY);
				DirectionX[size] = i;
				DirectionY[size] = num;
				SelIndex[size] = size;
				size++;
			}
		}
		//
		bool switched = false;
		do
		{
			switched = false;
			for(int i = 1; i < size; i++)
			{
				if(Sel[SelIndex[i - 1]] < Sel[SelIndex[i]])
				{
					int temp = SelIndex[i];
					SelIndex[i] = SelIndex[i - 1];
					SelIndex[i - 1] = temp;
					switched = true;
				}
			}
		}
		while(switched);
		//
		for(int i=0;i<size;i++)
		{
			if( DEBUG )
				cout << "sorted " << SelIndex[i] << " " << Sel[SelIndex[i]] << "--> " << DirectionX[SelIndex[i]] << " " << DirectionY[SelIndex[i]] << endl;
		}
		cout << ientry << " sorted " << SelIndex[0] << " " << Sel[SelIndex[0]] << "--> " << DirectionX[SelIndex[0]] << " " << DirectionY[SelIndex[0]] << " " << 7-DirectionX[SelIndex[0]]+2 << " " << DirectionY[SelIndex[0]]+2 << endl;
		
		bin = hpl2d->GetBin(7-DirectionX[SelIndex[0]]+2,DirectionY[SelIndex[0]]+2);
		hpl2d->AddBinContent(bin,1);
		/////
		PossibleCoincidence = TMath::LocMax(8,weightedBottom);
		PossibleCoincidence2 = TMath::LocMax(8,weightedTop);
		///
		if( weightedBottom[PossibleCoincidence] > weightedTop[PossibleCoincidence2] )
		{
			cout << ientry << " bottom signal arrives first " << 7-PossibleCoincidence+2 << " " << PossibleCoincidence2+2 << endl;
			binb = hplbot->GetBin(7-PossibleCoincidence+2,PossibleCoincidence2+2);
			hplbot->AddBinContent(binb,1);
		} 
		if( weightedTop[PossibleCoincidence2] > weightedBottom[PossibleCoincidence] )
		{
			cout << ientry << " Top signal arrives first " << 7-PossibleCoincidence+2 << " " << PossibleCoincidence2+2 << endl;
			bint = hpltop->GetBin(7-PossibleCoincidence+2,PossibleCoincidence2+2);
			hpltop->AddBinContent(bint,1);
		}
		if( weightedTop[PossibleCoincidence2] == weightedBottom[PossibleCoincidence] )
		{
			cout << ientry << " Both signal arrives first " << 7-PossibleCoincidence+2 << " " << PossibleCoincidence2+2 << endl;
			bineq = hpleq->GetBin(7-PossibleCoincidence+2,PossibleCoincidence2+2);
			hpleq->AddBinContent(bineq,1);
		}
		cout << ientry 
			 << " locate max. " 
			 << PossibleCoincidence 
			 << " " << weightedBottom[PossibleCoincidence] 
			 << " CH#" << 2*PossibleCoincidence 
			 << "-CH#" << 2*PossibleCoincidence+1 
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence+1)
			 << " locate max2. " 
			 << PossibleCoincidence2 
			 << " " << weightedTop[PossibleCoincidence2] 
			 << " CH#" << 2*PossibleCoincidence2+16
			 << "-CH#" << 2*PossibleCoincidence2+1+16
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence2+16)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence2+1+16)
			 << endl;
	}
	ps->NewPage();
	hpl2d->GetYaxis()->SetBinLabel(9,"Scin16");
	hpl2d->GetYaxis()->SetBinLabel(8,"Scin15");
	hpl2d->GetYaxis()->SetBinLabel(7,"Scin14");
	hpl2d->GetYaxis()->SetBinLabel(6,"Scin13");
	hpl2d->GetYaxis()->SetBinLabel(5,"Scin12");
	hpl2d->GetYaxis()->SetBinLabel(4,"Scin11");
	hpl2d->GetYaxis()->SetBinLabel(3,"Scin10");
	hpl2d->GetYaxis()->SetBinLabel(2,"Scin9");
	hpl2d->GetYaxis()->SetLabelSize(0.04);
	//
	hpl2d->GetXaxis()->SetBinLabel(9,"Scin1");
	hpl2d->GetXaxis()->SetBinLabel(8,"Scin2");
	hpl2d->GetXaxis()->SetBinLabel(7,"Scin3");
	hpl2d->GetXaxis()->SetBinLabel(6,"Scin4");
	hpl2d->GetXaxis()->SetBinLabel(5,"Scin5");
	hpl2d->GetXaxis()->SetBinLabel(4,"Scin6");
	hpl2d->GetXaxis()->SetBinLabel(3,"Scin7");
	hpl2d->GetXaxis()->SetBinLabel(2,"Scin8");
	hpl2d->GetXaxis()->SetLabelSize(0.04);
	hpl2d->SetStats(kFALSE);
	//
	//hpl2d->SetTitle(Form("Cosmic Event %i",ientry));
	hpl2d->GetXaxis()->SetTitleOffset(1.5);
	hpl2d->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
	hpl2d->GetYaxis()->SetTitleOffset(1.5);
	hpl2d->GetYaxis()->SetTitle("Top Layer (10 mm)");
	hpl2d->DrawCopy("col text");
	//
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hpl2d->GetZaxis()->SetRangeUser(min, max);
	//
	ps->NewPage();
	hpltop->GetYaxis()->SetBinLabel(9,"Scin16");
	hpltop->GetYaxis()->SetBinLabel(8,"Scin15");
	hpltop->GetYaxis()->SetBinLabel(7,"Scin14");
	hpltop->GetYaxis()->SetBinLabel(6,"Scin13");
	hpltop->GetYaxis()->SetBinLabel(5,"Scin12");
	hpltop->GetYaxis()->SetBinLabel(4,"Scin11");
	hpltop->GetYaxis()->SetBinLabel(3,"Scin10");
	hpltop->GetYaxis()->SetBinLabel(2,"Scin9");
	hpltop->GetYaxis()->SetLabelSize(0.04);
	//
	hpltop->GetXaxis()->SetBinLabel(9,"Scin1");
	hpltop->GetXaxis()->SetBinLabel(8,"Scin2");
	hpltop->GetXaxis()->SetBinLabel(7,"Scin3");
	hpltop->GetXaxis()->SetBinLabel(6,"Scin4");
	hpltop->GetXaxis()->SetBinLabel(5,"Scin5");
	hpltop->GetXaxis()->SetBinLabel(4,"Scin6");
	hpltop->GetXaxis()->SetBinLabel(3,"Scin7");
	hpltop->GetXaxis()->SetBinLabel(2,"Scin8");
	hpltop->GetXaxis()->SetLabelSize(0.04);
	hpltop->SetStats(kFALSE);
	hpltop->DrawCopy("col text");
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hpltop->GetZaxis()->SetRangeUser(min, max);
	//
	ps->NewPage();
	hplbot->GetYaxis()->SetBinLabel(9,"Scin16");
	hplbot->GetYaxis()->SetBinLabel(8,"Scin15");
	hplbot->GetYaxis()->SetBinLabel(7,"Scin14");
	hplbot->GetYaxis()->SetBinLabel(6,"Scin13");
	hplbot->GetYaxis()->SetBinLabel(5,"Scin12");
	hplbot->GetYaxis()->SetBinLabel(4,"Scin11");
	hplbot->GetYaxis()->SetBinLabel(3,"Scin10");
	hplbot->GetYaxis()->SetBinLabel(2,"Scin9");
	hplbot->GetYaxis()->SetLabelSize(0.04);
	//
	hplbot->GetXaxis()->SetBinLabel(9,"Scin1");
	hplbot->GetXaxis()->SetBinLabel(8,"Scin2");
	hplbot->GetXaxis()->SetBinLabel(7,"Scin3");
	hplbot->GetXaxis()->SetBinLabel(6,"Scin4");
	hplbot->GetXaxis()->SetBinLabel(5,"Scin5");
	hplbot->GetXaxis()->SetBinLabel(4,"Scin6");
	hplbot->GetXaxis()->SetBinLabel(3,"Scin7");
	hplbot->GetXaxis()->SetBinLabel(2,"Scin8");
	hplbot->GetXaxis()->SetLabelSize(0.04);
	hplbot->SetStats(kFALSE);
	hplbot->DrawCopy("col text");
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hplbot->GetZaxis()->SetRangeUser(min, max);
	//
	ps->NewPage();
	hpleq->GetYaxis()->SetBinLabel(9,"Scin16");
	hpleq->GetYaxis()->SetBinLabel(8,"Scin15");
	hpleq->GetYaxis()->SetBinLabel(7,"Scin14");
	hpleq->GetYaxis()->SetBinLabel(6,"Scin13");
	hpleq->GetYaxis()->SetBinLabel(5,"Scin12");
	hpleq->GetYaxis()->SetBinLabel(4,"Scin11");
	hpleq->GetYaxis()->SetBinLabel(3,"Scin10");
	hpleq->GetYaxis()->SetBinLabel(2,"Scin9");
	hpleq->GetYaxis()->SetLabelSize(0.04);
	//
	hpleq->GetXaxis()->SetBinLabel(9,"Scin1");
	hpleq->GetXaxis()->SetBinLabel(8,"Scin2");
	hpleq->GetXaxis()->SetBinLabel(7,"Scin3");
	hpleq->GetXaxis()->SetBinLabel(6,"Scin4");
	hpleq->GetXaxis()->SetBinLabel(5,"Scin5");
	hpleq->GetXaxis()->SetBinLabel(4,"Scin6");
	hpleq->GetXaxis()->SetBinLabel(3,"Scin7");
	hpleq->GetXaxis()->SetBinLabel(2,"Scin8");
	hpleq->GetXaxis()->SetLabelSize(0.04);
	hpleq->SetStats(kFALSE);
	hpleq->DrawCopy("col text");
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hpleq->GetZaxis()->SetRangeUser(min, max);
	
		for (Int_t biny=1;biny<=nbinsy;biny++) {
		for (Int_t binx=1;binx<=nbinsx;binx++) {
			ibin = hpl2d->GetBin(binx,biny);
			
			ibint = hpltop->GetBin(binx,biny);
			
			ibinb = hplbot->GetBin(binx,biny);
			
			ibineq = hpleq->GetBin(binx,biny);
			
			cout << "Bin[" << binx << "-" << biny << "]"
				 << " " << hpl2d->GetBinContent(ibin) 
				 << " " << hpltop->GetBinContent(ibint)
				 << " " << hplbot->GetBinContent(ibinb)
				 << " " << hpleq->GetBinContent(ibineq)
				 << endl;
		}
	}
	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ShowEventOnCRT_test(int ientry)
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
//	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
//	if (C) delete C;
	
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	
//	TCanvas *myCan = new TCanvas("Canvas","Show Event on CRT",200,10,700,500);
//	TPostScript *ps = new TPostScript(Form("CRT_Event_%i.ps",ientry),112);

	TH2I *hpl2d = new TH2I("histEv","Cosmic muon events",10,-1,9,10,-1,9);
	Int_t  nbinsx = hpl2d->GetXaxis()->GetNbins();
	Int_t  nbinsy = hpl2d->GetYaxis()->GetNbins();
	Int_t ibin,bin;

	for (Int_t biny=1;biny<=nbinsy;biny++) {
		for (Int_t binx=1;binx<=nbinsx;binx++) {
			ibin = hpl2d->GetBin(binx,biny);
			hpl2d->SetBinContent(ibin,0);
		}
	}

	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	int num, index = 0;
	Double_t ratio = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;
	data->GetMppc(ientry);
	Double_t Sel[64];
	Int_t SelIndex[64];
	Int_t DirectionX[64];
	Int_t DirectionY[64];
	Int_t size=0;
	//
	for(int i = 0; i < 8; i++)
	{
		for(int j = 8; j < 16; j++)
		{
			BarX = data->ADC(2*i) + data->ADC(2*i+1);
			BarY = data->ADC(2*j) + data->ADC(2*j+1);
			if(BarX==0||BarY==0) continue;
			//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
			num = abs(8-j);
			if(BarX==0||BarY==0) 
			{
				ratio = 0;
			} else{
				ratio = (Double_t) BarX/BarY;
			}
			if(ratio>1.)
			{
				ratio = (Double_t) BarY/BarX;
			}
			Sel[size] = ratio*(BarX+BarY);
			DirectionX[size] = 7 - i + 2;
			DirectionY[size] = num + 2;
			SelIndex[size] = size;
			size++;
		}
	}
	bool switched = false;
	do
	{
		switched = false;
		for(int i = 1; i < size; i++)
		{
			if(Sel[SelIndex[i - 1]] < Sel[SelIndex[i]])
			{
				int temp = SelIndex[i];
				SelIndex[i] = SelIndex[i - 1];
				SelIndex[i - 1] = temp;
				switched = true;
			}
		}
	}
	while(switched);
	//
	//for(int i=0;i<size;i++)
	//{
	//		cout << "sorted " << SelIndex[i] << " " << Sel[SelIndex[i]] << "--> " << DirectionX[SelIndex[i]] << " " << DirectionY[SelIndex[i]] << endl;
	//}
	bin = hpl2d->GetBin(DirectionX[SelIndex[0]],DirectionY[SelIndex[0]]);
	hpl2d->AddBinContent(bin,1);
	
	TCanvas *c = new TCanvas("Canvas","Show Event on CRT",200,10,700,500);
	c->cd();
	hpl2d->GetYaxis()->SetBinLabel(9,"Scin16");
	hpl2d->GetYaxis()->SetBinLabel(8,"Scin15");
	hpl2d->GetYaxis()->SetBinLabel(7,"Scin14");
	hpl2d->GetYaxis()->SetBinLabel(6,"Scin13");
	hpl2d->GetYaxis()->SetBinLabel(5,"Scin12");
	hpl2d->GetYaxis()->SetBinLabel(4,"Scin11");
	hpl2d->GetYaxis()->SetBinLabel(3,"Scin10");
	hpl2d->GetYaxis()->SetBinLabel(2,"Scin9");
	hpl2d->GetYaxis()->SetLabelSize(0.04);
	//
	hpl2d->GetXaxis()->SetBinLabel(9,"Scin1");
	hpl2d->GetXaxis()->SetBinLabel(8,"Scin2");
	hpl2d->GetXaxis()->SetBinLabel(7,"Scin3");
	hpl2d->GetXaxis()->SetBinLabel(6,"Scin4");
	hpl2d->GetXaxis()->SetBinLabel(5,"Scin5");
	hpl2d->GetXaxis()->SetBinLabel(4,"Scin6");
	hpl2d->GetXaxis()->SetBinLabel(3,"Scin7");
	hpl2d->GetXaxis()->SetBinLabel(2,"Scin8");
	hpl2d->GetXaxis()->SetLabelSize(0.04);
	hpl2d->SetStats(kFALSE);
	
	hpl2d->SetTitle(Form("Cosmic Event %i",ientry));
	hpl2d->GetXaxis()->SetTitleOffset(1.5);
	hpl2d->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
	hpl2d->GetYaxis()->SetTitleOffset(1.5);
	hpl2d->GetYaxis()->SetTitle("Top Layer (10 mm)");
	hpl2d->DrawCopy("col text");
	c->Update();
	c->Modified();
	c->Update();
	hpl2d->GetZaxis()->SetRangeUser(min, max);
//	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ShowSingleEventOnCRT_Fraction(int ientry)
{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	TH2I *mt = (TH2I*) gROOT->FindObject("histEvTop");
	if(mt) delete mt;
	TH2I *mb = (TH2I*) gROOT->FindObject("histEvBot");
	if(mb) delete mb;
	

	TCanvas *myCan = new TCanvas("Canvas","Show Event on CRT",200,10,700,500);
	TPostScript *ps = new TPostScript(Form("CRT_SingleEvent_Fraction_%i.ps",ientry),112);

	TH2I *hpl2d = new TH2I("histEv","Cosmic muon event",10,-1,9,10,-1,9);
	TH2I *hpltop = new TH2I("histTop","Cosmic muon event (Top Layer Signal Arrives First)",10,-1,9,10,-1,9);
	TH2I *hplbot = new TH2I("histBot","Cosmic muon event (Bottom Layer Signal Arrives First)",10,-1,9,10,-1,9);
	TH2I *hpleq = new TH2I("histEq","Cosmic muon event (Both Signal has same strength)",10,-1,9,10,-1,9);

	Int_t  nbinsx = hpl2d->GetXaxis()->GetNbins();
	Int_t  nbinsy = hpl2d->GetYaxis()->GetNbins();
	Int_t ibin,bin,ibint,bint,ibinb,binb,ibineq,bineq;
	for (Int_t biny=1;biny<=nbinsy;biny++) {
		for (Int_t binx=1;binx<=nbinsx;binx++) {
			ibin = hpl2d->GetBin(binx,biny);
			hpl2d->SetBinContent(ibin,0);
			ibint = hpltop->GetBin(binx,biny);
			hpltop->SetBinContent(ibint,0);
			ibinb = hplbot->GetBin(binx,biny);
			hplbot->SetBinContent(ibinb,0);
			ibineq = hpleq->GetBin(binx,biny);
			hpleq->SetBinContent(ibineq,0);
		}
	}
	//
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	int num;
	Double_t ratio = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;
	std::vector<int> v_InCoincidence;
	std::vector<double> sel;
	TBox *Box[10];
	UInt_t PossibleCoincidence;
	UInt_t PossibleCoincidence2;
	Double_t ratioBottom[8];
	Double_t weightedBottom[8];
	Double_t ratioTop[8];
	Double_t weightedTop[8];
	//
	Double_t Sel[64];
	Int_t SelIndex[64];
	Int_t DirectionX[64];
	Int_t DirectionY[64];
	Int_t size;
//	for( Long64_t ientry=0; ientry<entries; ++ientry )
//	{
		data->GetMppc(ientry);
//		if(ientry>1000) break;
		size = 0;
		for(int i = 0; i < 8; i++)
		{
			if (data->ADC(2*i) == 0 || data->ADC(2*i+1) ==0)
			{
				ratioBottom[i]    = 0;
				weightedBottom[i] = 0;
			}else{
				ratioBottom[i] = (Double_t)data->ADC(2*i)/data->ADC(2*i+1);
				if(ratioBottom[i]>1.) 
				{
					ratioBottom[i] = (Double_t)data->ADC(2*i+1)/data->ADC(2*i);
				}
				weightedBottom[i] = ratioBottom[i]*(data->ADC(2*i)+data->ADC(2*i+1));
			}
			if (data->ADC(2*i+16) == 0 || data->ADC(2*i+1+16) ==0)
			{
				ratioTop[i]    = 0;
				weightedTop[i] = 0;
			}else{
				ratioTop[i] = (Double_t)data->ADC(2*i+16)/data->ADC(2*i+1+16);
				if(ratioTop[i]>1.) 
				{
					ratioTop[i] = (Double_t)data->ADC(2*i+1+16)/data->ADC(2*i+16);
				}
				weightedTop[i] = ratioTop[i]*(data->ADC(2*i+1+16)+data->ADC(2*i+1+16));
			}
			///////
			for(int j = 8; j < 16; j++)
			{
				BarX = data->ADC(2*i) + data->ADC(2*i+1);
				BarY = data->ADC(2*j) + data->ADC(2*j+1);
				if(BarX==0||BarY==0) continue;
				//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
				num = abs(8-j);
				ratio = (Double_t) BarX/BarY;
				if(ratio>1.)
				{
					ratio = (Double_t) BarY/BarX;
				}
				if( DEBUG ) 
					cout << "---> " << ratio*(BarX+BarY) << " " << i << " " << num << " " << size << " :: " << 8*i+num<< endl;
				Sel[size] = ratio*(BarX+BarY);
				DirectionX[size] = i;
				DirectionY[size] = num;
				SelIndex[size] = size;
				size++;
			}
		}
		//
		bool switched = false;
		do
		{
			switched = false;
			for(int i = 1; i < size; i++)
			{
				if(Sel[SelIndex[i - 1]] < Sel[SelIndex[i]])
				{
					int temp = SelIndex[i];
					SelIndex[i] = SelIndex[i - 1];
					SelIndex[i - 1] = temp;
					switched = true;
				}
			}
		}
		while(switched);
		//
		for(int i=0;i<size;i++)
		{
			if( DEBUG )
				cout << "sorted " << SelIndex[i] << " " << Sel[SelIndex[i]] << "--> " << DirectionX[SelIndex[i]] << " " << DirectionY[SelIndex[i]] << endl;
		}
		cout << ientry << " sorted " << SelIndex[0] << " " << Sel[SelIndex[0]] << "--> " << DirectionX[SelIndex[0]] << " " << DirectionY[SelIndex[0]] << " " << 7-DirectionX[SelIndex[0]]+2 << " " << DirectionY[SelIndex[0]]+2 << endl;
		
		bin = hpl2d->GetBin(7-DirectionX[SelIndex[0]]+2,DirectionY[SelIndex[0]]+2);
		hpl2d->AddBinContent(bin,1);
		/////
		PossibleCoincidence = TMath::LocMax(8,weightedBottom);
		PossibleCoincidence2 = TMath::LocMax(8,weightedTop);
		///
		if( weightedBottom[PossibleCoincidence] > weightedTop[PossibleCoincidence2] )
		{
			cout << ientry << " bottom signal arrives first " << 7-PossibleCoincidence+2 << " " << PossibleCoincidence2+2 << endl;
			binb = hplbot->GetBin(7-PossibleCoincidence+2,PossibleCoincidence2+2);
			hplbot->AddBinContent(binb,1);
		} 
		if( weightedTop[PossibleCoincidence2] > weightedBottom[PossibleCoincidence] )
		{
			cout << ientry << " Top signal arrives first " << 7-PossibleCoincidence+2 << " " << PossibleCoincidence2+2 << endl;
			bint = hpltop->GetBin(7-PossibleCoincidence+2,PossibleCoincidence2+2);
			hpltop->AddBinContent(bint,1);
		}
		if( weightedTop[PossibleCoincidence2] == weightedBottom[PossibleCoincidence] )
		{
			cout << ientry << " Both signal arrives first " << 7-PossibleCoincidence+2 << " " << PossibleCoincidence2+2 << endl;
			bineq = hpleq->GetBin(7-PossibleCoincidence+2,PossibleCoincidence2+2);
			hpleq->AddBinContent(bineq,1);
		}
		cout << ientry 
			 << " locate max. " 
			 << PossibleCoincidence 
			 << " " << weightedBottom[PossibleCoincidence] 
			 << " CH#" << 2*PossibleCoincidence 
			 << "-CH#" << 2*PossibleCoincidence+1 
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence+1)
			 << " locate max2. " 
			 << PossibleCoincidence2 
			 << " " << weightedTop[PossibleCoincidence2] 
			 << " CH#" << 2*PossibleCoincidence2+16
			 << "-CH#" << 2*PossibleCoincidence2+1+16
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence2+16)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence2+1+16)
			 << endl;
//	}
	ps->NewPage();
	hpl2d->GetYaxis()->SetBinLabel(9,"Scin16");
	hpl2d->GetYaxis()->SetBinLabel(8,"Scin15");
	hpl2d->GetYaxis()->SetBinLabel(7,"Scin14");
	hpl2d->GetYaxis()->SetBinLabel(6,"Scin13");
	hpl2d->GetYaxis()->SetBinLabel(5,"Scin12");
	hpl2d->GetYaxis()->SetBinLabel(4,"Scin11");
	hpl2d->GetYaxis()->SetBinLabel(3,"Scin10");
	hpl2d->GetYaxis()->SetBinLabel(2,"Scin9");
	hpl2d->GetYaxis()->SetLabelSize(0.04);
	//
	hpl2d->GetXaxis()->SetBinLabel(9,"Scin1");
	hpl2d->GetXaxis()->SetBinLabel(8,"Scin2");
	hpl2d->GetXaxis()->SetBinLabel(7,"Scin3");
	hpl2d->GetXaxis()->SetBinLabel(6,"Scin4");
	hpl2d->GetXaxis()->SetBinLabel(5,"Scin5");
	hpl2d->GetXaxis()->SetBinLabel(4,"Scin6");
	hpl2d->GetXaxis()->SetBinLabel(3,"Scin7");
	hpl2d->GetXaxis()->SetBinLabel(2,"Scin8");
	hpl2d->GetXaxis()->SetLabelSize(0.04);
	hpl2d->SetStats(kFALSE);
	//
	//hpl2d->SetTitle(Form("Cosmic Event %i",ientry));
	hpl2d->GetXaxis()->SetTitleOffset(1.5);
	hpl2d->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
	hpl2d->GetYaxis()->SetTitleOffset(1.5);
	hpl2d->GetYaxis()->SetTitle("Top Layer (10 mm)");
	hpl2d->DrawCopy("col text");
	//
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hpl2d->GetZaxis()->SetRangeUser(min, max);
	//
	ps->NewPage();
	hpltop->GetYaxis()->SetBinLabel(9,"Scin16");
	hpltop->GetYaxis()->SetBinLabel(8,"Scin15");
	hpltop->GetYaxis()->SetBinLabel(7,"Scin14");
	hpltop->GetYaxis()->SetBinLabel(6,"Scin13");
	hpltop->GetYaxis()->SetBinLabel(5,"Scin12");
	hpltop->GetYaxis()->SetBinLabel(4,"Scin11");
	hpltop->GetYaxis()->SetBinLabel(3,"Scin10");
	hpltop->GetYaxis()->SetBinLabel(2,"Scin9");
	hpltop->GetYaxis()->SetLabelSize(0.04);
	//
	hpltop->GetXaxis()->SetBinLabel(9,"Scin1");
	hpltop->GetXaxis()->SetBinLabel(8,"Scin2");
	hpltop->GetXaxis()->SetBinLabel(7,"Scin3");
	hpltop->GetXaxis()->SetBinLabel(6,"Scin4");
	hpltop->GetXaxis()->SetBinLabel(5,"Scin5");
	hpltop->GetXaxis()->SetBinLabel(4,"Scin6");
	hpltop->GetXaxis()->SetBinLabel(3,"Scin7");
	hpltop->GetXaxis()->SetBinLabel(2,"Scin8");
	hpltop->GetXaxis()->SetLabelSize(0.04);
	hpltop->SetStats(kFALSE);
	hpltop->DrawCopy("col text");
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hpltop->GetZaxis()->SetRangeUser(min, max);
	//
	ps->NewPage();
	hplbot->GetYaxis()->SetBinLabel(9,"Scin16");
	hplbot->GetYaxis()->SetBinLabel(8,"Scin15");
	hplbot->GetYaxis()->SetBinLabel(7,"Scin14");
	hplbot->GetYaxis()->SetBinLabel(6,"Scin13");
	hplbot->GetYaxis()->SetBinLabel(5,"Scin12");
	hplbot->GetYaxis()->SetBinLabel(4,"Scin11");
	hplbot->GetYaxis()->SetBinLabel(3,"Scin10");
	hplbot->GetYaxis()->SetBinLabel(2,"Scin9");
	hplbot->GetYaxis()->SetLabelSize(0.04);
	//
	hplbot->GetXaxis()->SetBinLabel(9,"Scin1");
	hplbot->GetXaxis()->SetBinLabel(8,"Scin2");
	hplbot->GetXaxis()->SetBinLabel(7,"Scin3");
	hplbot->GetXaxis()->SetBinLabel(6,"Scin4");
	hplbot->GetXaxis()->SetBinLabel(5,"Scin5");
	hplbot->GetXaxis()->SetBinLabel(4,"Scin6");
	hplbot->GetXaxis()->SetBinLabel(3,"Scin7");
	hplbot->GetXaxis()->SetBinLabel(2,"Scin8");
	hplbot->GetXaxis()->SetLabelSize(0.04);
	hplbot->SetStats(kFALSE);
	hplbot->DrawCopy("col text");
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hplbot->GetZaxis()->SetRangeUser(min, max);
	//
	ps->NewPage();
	hpleq->GetYaxis()->SetBinLabel(9,"Scin16");
	hpleq->GetYaxis()->SetBinLabel(8,"Scin15");
	hpleq->GetYaxis()->SetBinLabel(7,"Scin14");
	hpleq->GetYaxis()->SetBinLabel(6,"Scin13");
	hpleq->GetYaxis()->SetBinLabel(5,"Scin12");
	hpleq->GetYaxis()->SetBinLabel(4,"Scin11");
	hpleq->GetYaxis()->SetBinLabel(3,"Scin10");
	hpleq->GetYaxis()->SetBinLabel(2,"Scin9");
	hpleq->GetYaxis()->SetLabelSize(0.04);
	//
	hpleq->GetXaxis()->SetBinLabel(9,"Scin1");
	hpleq->GetXaxis()->SetBinLabel(8,"Scin2");
	hpleq->GetXaxis()->SetBinLabel(7,"Scin3");
	hpleq->GetXaxis()->SetBinLabel(6,"Scin4");
	hpleq->GetXaxis()->SetBinLabel(5,"Scin5");
	hpleq->GetXaxis()->SetBinLabel(4,"Scin6");
	hpleq->GetXaxis()->SetBinLabel(3,"Scin7");
	hpleq->GetXaxis()->SetBinLabel(2,"Scin8");
	hpleq->GetXaxis()->SetLabelSize(0.04);
	hpleq->SetStats(kFALSE);
	hpleq->DrawCopy("col text");
	gPad->Modified();
	gPad->Update();
	myCan->Modified();
	myCan->Update();
	hpleq->GetZaxis()->SetRangeUser(min, max);
	
		for (Int_t biny=1;biny<=nbinsy;biny++) {
		for (Int_t binx=1;binx<=nbinsx;binx++) {
			ibin = hpl2d->GetBin(binx,biny);
			
			ibint = hpltop->GetBin(binx,biny);
			
			ibinb = hplbot->GetBin(binx,biny);
			
			ibineq = hpleq->GetBin(binx,biny);
			
			cout << "Bin[" << binx << "-" << biny << "]"
				 << " " << hpl2d->GetBinContent(ibin) 
				 << " " << hpltop->GetBinContent(ibint)
				 << " " << hplbot->GetBinContent(ibinb)
				 << " " << hpleq->GetBinContent(ibineq)
				 << endl;
		}
	}

	
	
	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void GetShowerOnCRT()
{
	cout << "Here we are getting list of shower events!!!!" << endl;
	TH2I *m = (TH2I*) gROOT->FindObject("histEv");
	if(m) delete m;
	

	TH2I *hpl2d = new TH2I("histEv","Cosmic muon event",10,-1,9,10,-1,9);
	Int_t  nbinsx = hpl2d->GetXaxis()->GetNbins();
	Int_t  nbinsy = hpl2d->GetYaxis()->GetNbins();
	Int_t ibin,bin;

	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	UInt_t BarX, BarY;
	UInt_t adc_value;
	int count1, count2, count3;
	int num, index = 0;
	Double_t ratio = 0.;
	Double_t ratio_t = 0.;
	Double_t ratio_b = 0.;
	const Double_t min = 0.;
	const Double_t max = 1.1;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		count1 = count2 = count3 = 0;
		data->GetMppc(ientry);
		//
		for(int i = 0; i < 8; i++)
		{
			for(int j = 8; j < 16; j++)
			{
				index++;
				BarX = data->ADC(2*i) + data->ADC(2*i+1);
				BarY = data->ADC(2*j) + data->ADC(2*j+1);
				//cout << 2*i << " " << 2*i+1 << " " << 2*j << " " << 2*j+1 << " " << num << endl;
				num = abs(8-j);
				if(BarX==0||BarY==0) 
				{
					ratio = 0;
				} else{
					ratio = (Double_t) BarX/BarY;
				}
				if(ratio>1.)
				{
					ratio = (Double_t) BarY/BarX;
				}
				//hpl2d->SetBinContent(i+2,num+2,ratio*(BarX+BarY));
				hpl2d->SetBinContent(7-i+2,num+2,ratio*(BarX+BarY));
			}
		}
		for (Int_t biny=1;biny<=nbinsy;biny++)
		{
			for (Int_t binx=1;binx<=nbinsx;binx++)
			{
				ibin = hpl2d->GetBin(binx,biny);
				adc_value = hpl2d->GetBinContent(ibin);
				if(adc_value>1000) count1++;
				if(adc_value>1500) count2++;
				if(adc_value>2000) count3++;
			}
		}
		if(count1 > 1) cout << "Shower!!! " << ientry << " --> hits " << count1 << " " << count2 << " " << count3 << endl;
	}	
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void CheckCoincidenceX(){
	TPostScript *ps = new TPostScript("Coincidence.ps",112);
	TH1F *histos[100];
	TH1F *h_CAENCoin[100];
	for(int i = 0; i < 100; i++) 
	{
		histos[i] = new TH1F(Form("histos[%i]",i),Form("Event[%i];Channel Numner;ADC",i),32,-0.5,31.5);
		h_CAENCoin[i] = new TH1F(Form("h_CAENCoin[%i]",i),Form("Event[%i];Channel Numner;ADC",i),32,-0.5,31.5);
	}
	TCanvas *myCan = new TCanvas();
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	cout << "Number of entries " << entries << endl;
	UInt_t BarX, BarY;
	Int_t num;
	Int_t InCoincidence;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
		//
		if(ientry<100){
			for (int i = 0; i < 32; i++ )
			{
				cout << i << " " << data->ADC(i) << " " << data->Coincidence() << " " << endl;
				InCoincidence = ConvertToBinary(ientry,data->Coincidence(),16);
				cout << endl;
				histos[ientry]->SetBinContent(i+1,data->ADC(i));
			}
			if(DEBUG) cout << "--> " << InCoincidence << endl;
			if(InCoincidence % 2 == 0)
			{
				h_CAENCoin[ientry]->SetBinContent(InCoincidence,data->ADC(InCoincidence-1));
				h_CAENCoin[ientry]->SetBinContent(InCoincidence+1,data->ADC(InCoincidence));
			}else
			{
				h_CAENCoin[ientry]->SetBinContent(InCoincidence,data->ADC(InCoincidence-1));
				h_CAENCoin[ientry]->SetBinContent(InCoincidence-1,data->ADC(InCoincidence-2));
			}
		}
	}
	for(int i = 0; i < 100; i++)
	{
		histos[i]->Draw();
		h_CAENCoin[i]->SetFillColor(kRed);
		h_CAENCoin[i]->SetFillStyle(3010);
		h_CAENCoin[i]->Draw("same");
		myCan->Update();
		ps->NewPage();
		//histos[i]->SetFillColor(kRed);
		//histos[i]->SetFillStyle(3010);
	}
	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void CheckCoincidence(){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TPostScript *ps = new TPostScript("Coincidence.ps",112);
	TH1F *histos[100];
	TH1F *h_CAENCoin[100];
	TH1F *h_PosCoinBottom[100];
	TH1F *h_PosCoinTop[100];
	Int_t ScintillatorX[100];
	Int_t ScintillatorY[100];
	
	TH2F *hCoinPlane = new TH2F("hCoinPlane","Cosmic Muon Event in Coincidence;15 mm Scintillator Bar;10 mm Scintillator Bar",20,0,8,20,0,8);

	for(int i = 0; i < 100; i++) 
	{
		histos[i] = new TH1F(Form("histos[%i]",i),Form("Cosmic Muon Event %i;Channel Number;ADC",i),32,-0.5,31.5);
		h_CAENCoin[i] = new TH1F(Form("h_CAENCoin[%i]",i),Form("Cosmic Muon Event %i;Channel Number;ADC",i),32,-0.5,31.5);
		h_PosCoinBottom[i] = new TH1F(Form("h_PosCoinBottom[%i]",i),Form("Cosmic Muon Event %i;Channel Number;ADC",i),32,-0.5,31.5);
		h_PosCoinTop[i] = new TH1F(Form("h_PosCoinTop[%i]",i),Form("Cosmic Muon Event %i;Channel Number;ADC",i),32,-0.5,31.5);
	}
	
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Check Coincidence",200,10,700,500);
	
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	cout << "Number of entries " << entries << endl;
	UInt_t BarX, BarY;
	Int_t num;
	Int_t InCoincidence;
	std::vector<int> v_InCoincidence;
	UInt_t PossibleCoincidence;
	UInt_t PossibleCoincidence2;
	UInt_t MaxAdcInPosCoincidence;
	UInt_t MinAdcInPosCoincidence;
	UInt_t MaxAdcInPosCoincidence2;
	UInt_t MinAdcInPosCoincidence2;

	Double_t ratioBottom[8];
	Double_t weightedBottom[8];
	Double_t ratioTop[8];
	Double_t weightedTop[8];
	TBox *LayerX, *LayerY;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
		if(ientry<100){
			for (int i=0;i<8;i++)
			{
				//cout << "even " << 2*i << " " << event1_adc[2*i] << endl;
				//cout << "odd " << 2*i+1 << " "<< event1_adc[2*i+1] << endl;
				if (data->ADC(2*i) == 0 || data->ADC(2*i+1) ==0)
				{
					ratioBottom[i]    = 0;
					weightedBottom[i] = 0;
				}else{
					ratioBottom[i] = (Double_t)data->ADC(2*i)/data->ADC(2*i+1);
					weightedBottom[i] = ratioBottom[i]*data->ADC(2*i+1);
					if(ratioBottom[i]>1.) {
						ratioBottom[i] = (Double_t)data->ADC(2*i+1)/data->ADC(2*i);
						weightedBottom[i] = ratioBottom[i]*data->ADC(2*i);
					}
				}
				if(DEBUG) 
					cout << "Bottom ----> ratio for CH#" << 2*i << "-CH#" << 2*i+1 << " " << ratioBottom[i] << " " << weightedBottom[i] << " " << endl;
				//
				if (data->ADC(2*i+16) == 0 || data->ADC(2*i+1+16) ==0)
				{
					ratioTop[i]    = 0;
					weightedTop[i] = 0;
				}else{
					ratioTop[i] = (Double_t)data->ADC(2*i+16)/data->ADC(2*i+1+16);
					weightedTop[i] = ratioTop[i]*data->ADC(2*i+1+16);
					if(ratioTop[i]>1.) {
						ratioTop[i] = (Double_t)data->ADC(2*i+1+16)/data->ADC(2*i+16);
						weightedTop[i] = ratioTop[i]*data->ADC(2*i+16);
					}
				}
				if(DEBUG) 
					cout << "Top ----> ratio for CH#" << 2*i+16 << "-CH#" << 2*i+1+16 << " " << ratioTop[i] << " " << weightedTop[i] << " " << endl;
			}
			PossibleCoincidence = TMath::LocMax(8,weightedBottom);
			MaxAdcInPosCoincidence = TMath::Max(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
			MinAdcInPosCoincidence = TMath::Min(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
			if(DEBUG) 
			cout << ientry 
				 << " locate max. " 
				 << PossibleCoincidence 
				 << " " << weightedBottom[PossibleCoincidence] 
				 << " CH#" << 2*PossibleCoincidence 
				 << "-CH#" << 2*PossibleCoincidence+1 
				 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence)
				 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence+1)
				 << " maxADC " << MaxAdcInPosCoincidence
				 << " minADC " << MinAdcInPosCoincidence
				 << endl;
			
			ScintillatorY[ientry] = PossibleCoincidence;
			
			Long64_t TwoChannelsInPosCoincidence_Plane1 = (2*PossibleCoincidence+2*PossibleCoincidence+1)/2;
			if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane1 " << TwoChannelsInPosCoincidence_Plane1 << endl;
			//
			PossibleCoincidence2 = TMath::LocMax(8,weightedTop);
			MaxAdcInPosCoincidence2 = TMath::Max(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
			MinAdcInPosCoincidence2 = TMath::Min(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
			if(DEBUG) 
			cout << ientry 
				 << " locate max2. " 
				 << PossibleCoincidence2 
				 << " " << weightedTop[PossibleCoincidence2] 
				 << " CH#" << 2*PossibleCoincidence2+16
				 << "-CH#" << 2*PossibleCoincidence2+1+16
				 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence2+16)
				 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence2+1+16)
				 << " maxADC " << MaxAdcInPosCoincidence2
				 << " minADC " << MinAdcInPosCoincidence2
				 << endl;


			hCoinPlane->Fill(PossibleCoincidence2+1,PossibleCoincidence+1);
			
			ScintillatorX[ientry] = PossibleCoincidence2;

			Long64_t TwoChannelsInPosCoincidence_Plane2 = (2*PossibleCoincidence2+16+2*PossibleCoincidence2+1+16)/2;
			if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane2 " << TwoChannelsInPosCoincidence_Plane2 << endl;

			for (int i = 0; i < 32; i++ )
			{
				cout << i << " " << data->ADC(i) << " " << data->Coincidence() << " " << endl;
				//if(i==31) InCoincidence = ConvertToBinary(ientry,data->Coincidence(),16);
				if(i==31) v_InCoincidence = ConvertToBinary2nd(ientry,data->Coincidence(),16);
				cout << endl;
				histos[ientry]->SetBinContent(i+1,data->ADC(i));
			}
			//if(DEBUG) cout << "Coincidence " << InCoincidence << " " << PossibleCoincidence << " " << PossibleCoincidence2 << endl;
			//if(DEBUG) cout << "CAEN coincidence pl " << InCoincidence << " " << 2*InCoincidence << " " << 2*InCoincidence+1 << endl;
			for(int i=0; i<v_InCoincidence.size();i++)
			{
				if(DEBUG) cout << "Coincidence " << v_InCoincidence[i] << " " << PossibleCoincidence << " " << PossibleCoincidence2 << endl;
				if(DEBUG) cout << "CAEN coincidence pl " << v_InCoincidence[i] << " " << 2*v_InCoincidence[i] << " " << 2*v_InCoincidence[i]+1 << endl;
				h_CAENCoin[ientry]->SetBinContent(2*v_InCoincidence[i]-1,data->ADC(2*v_InCoincidence[i]-2));
				h_CAENCoin[ientry]->SetBinContent(2*v_InCoincidence[i],data->ADC(2*v_InCoincidence[i]-1));
			}
			h_PosCoinBottom[ientry]->SetBinContent(TwoChannelsInPosCoincidence_Plane1+1,data->ADC(TwoChannelsInPosCoincidence_Plane1));
			h_PosCoinBottom[ientry]->SetBinContent(TwoChannelsInPosCoincidence_Plane1+2,data->ADC(TwoChannelsInPosCoincidence_Plane1+1));
			h_PosCoinTop[ientry]->SetBinContent(TwoChannelsInPosCoincidence_Plane2+1,data->ADC(TwoChannelsInPosCoincidence_Plane2));
			h_PosCoinTop[ientry]->SetBinContent(TwoChannelsInPosCoincidence_Plane2+2,data->ADC(TwoChannelsInPosCoincidence_Plane2+1));
		}// loop entry<100
	}
	for(int i = 0; i < 100; i++)
	{
		histos[i]->SetMaximum(5000);
		histos[i]->Draw();
		histos[i]->GetYaxis()->SetTitleOffset(1.3);
		h_CAENCoin[i]->SetFillColor(kRed);
		h_CAENCoin[i]->SetFillStyle(1001);
		h_CAENCoin[i]->Draw("same");
		h_PosCoinBottom[i]->SetFillColor(kGreen);
		h_PosCoinBottom[i]->SetFillStyle(3012);
		h_PosCoinBottom[i]->Draw("same");
		h_PosCoinTop[i]->SetFillColor(kGreen);
		h_PosCoinTop[i]->SetFillStyle(3022);
		h_PosCoinTop[i]->Draw("same");
		TLine *line = new TLine(15.5,0,15.5,5000);
		line->SetLineColor(kBlue);
		line->SetLineWidth(2);
		line->Draw("same");
		TText *t = new TText(7.5,4500,"Bottom Layer (15mm Bars)");
		t->SetTextAlign(22);
		t->SetTextColor(kBlue);
		t->SetTextFont(43);
		t->SetTextSize(14);
		t->Draw("same");
		TText *t2 = new TText(23,4500,"Top Layer (10mm Bars)");
		t2->SetTextAlign(22);
		t2->SetTextColor(kBlue);
		t2->SetTextFont(43);
		t2->SetTextSize(14);
		t2->Draw("same");
		myCan->Update();
		ps->NewPage();
		//myCan->Clear();
		TH1 *frame = myCan->DrawFrame(0,0,1850,1850);
		frame->GetXaxis()->SetBinLabel(80,"Scin8");
		frame->GetXaxis()->SetBinLabel(200,"Scin7");
		frame->GetXaxis()->SetBinLabel(320,"Scin6");
		frame->GetXaxis()->SetBinLabel(440,"Scin5");
		frame->GetXaxis()->SetBinLabel(560,"Scin4");
		frame->GetXaxis()->SetBinLabel(680,"Scin3");
		frame->GetXaxis()->SetBinLabel(800,"Scin2");
		frame->GetXaxis()->SetBinLabel(920,"Scin1");
		frame->GetXaxis()->SetLabelSize(0.04);

		frame->GetYaxis()->Set(8,frame->GetMinimum(),frame->GetMaximum());
		frame->GetYaxis()->SetBinLabel(8,"Scin16");
		frame->GetYaxis()->SetBinLabel(7,"Scin15");
		frame->GetYaxis()->SetBinLabel(6,"Scin14");
		frame->GetYaxis()->SetBinLabel(5,"Scin13");
		frame->GetYaxis()->SetBinLabel(4,"Scin12");
		frame->GetYaxis()->SetBinLabel(3,"Scin11");
		frame->GetYaxis()->SetBinLabel(2,"Scin10");
		frame->GetYaxis()->SetBinLabel(1,"Scin9");
		frame->GetYaxis()->SetLabelSize(0.04);

		frame->GetXaxis()->SetTitleOffset(1.6);
		frame->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
		frame->GetYaxis()->SetTitleOffset(1.6);
		frame->GetYaxis()->SetTitle("Top Layer (10 mm)");
		frame->SetTitle("CRT Module (Scintillators in Coincidence)");
		myCan->Update();
		//
		LayerY = new TBox(1850-230*ScintillatorY[i],0,1850-230*(ScintillatorY[i]+1),1850);
		LayerY->SetFillColor(kGray+1);
		LayerY->SetFillStyle(3006);
		LayerY->SetLineColor(kBlue);
		LayerY->SetLineWidth(2);
		LayerY->Draw();
		LayerX = new TBox(0,230*ScintillatorX[i],1850,230*ScintillatorX[i]+230);
		LayerX->SetFillColor(kGray+1);
		LayerX->SetFillStyle(3007);
		LayerX->SetLineColor(kRed);
		LayerX->SetLineWidth(2);
		LayerX->Draw("same");
		myCan->Update();
		ps->NewPage();
		//histos[i]->SetFillColor(kRed);
		//histos[i]->SetFillStyle(3010);
	}
		hCoinPlane->Draw();
		myCan->Update();
		ps->NewPage();
	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void CheckCoincidenceAll(){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TPostScript *ps = new TPostScript("Coincidence.ps",112);
	
	TH2F *hCoinPlane = new TH2F("hCoinPlane","Cosmic Muon Event in Coincidence;15 mm Scintillator Bar;10 mm Scintillator Bar",20,0,8,20,0,8);

	
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Check Coincidence All",200,10,700,500);

	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	cout << "Number of entries " << entries << endl;
	UInt_t BarX, BarY;
	Int_t num;
	Int_t InCoincidence;
	std::vector<int> v_InCoincidence;
	UInt_t PossibleCoincidence;
	UInt_t PossibleCoincidence2;
	UInt_t MaxAdcInPosCoincidence;
	UInt_t MinAdcInPosCoincidence;
	UInt_t MaxAdcInPosCoincidence2;
	UInt_t MinAdcInPosCoincidence2;
	
	Int_t total = 0, multip = 0;
	
	Double_t ratioBottom[8];
	Double_t weightedBottom[8];
	Double_t ratioTop[8];
	Double_t weightedTop[8];
	TBox *LayerX, *LayerY;
	for( Int_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
		TH1F *histos = new TH1F(Form("histos[%i]",ientry),Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);
		TH1F *h_CAENCoin = new TH1F(Form("h_CAENCoin[%i]",ientry),Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);
		TH1F *h_PosCoinBottom = new TH1F(Form("h_PosCoinBottom[%i]",ientry),Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);
		TH1F *h_PosCoinTop = new TH1F(Form("h_PosCoinTop[%i]",ientry),Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);
		total++;
		for (int i=0;i<8;i++)
		{
			//cout << "even " << 2*i << " " << event1_adc[2*i] << endl;
			//cout << "odd " << 2*i+1 << " "<< event1_adc[2*i+1] << endl;
			if (data->ADC(2*i) == 0 || data->ADC(2*i+1) ==0)
			{
				ratioBottom[i]    = 0;
				weightedBottom[i] = 0;
			}else{
				ratioBottom[i] = (Double_t)data->ADC(2*i)/data->ADC(2*i+1);
				weightedBottom[i] = ratioBottom[i]*data->ADC(2*i+1);
				if(ratioBottom[i]>1.) {
					ratioBottom[i] = (Double_t)data->ADC(2*i+1)/data->ADC(2*i);
					weightedBottom[i] = ratioBottom[i]*data->ADC(2*i);
				}
			}
			if(DEBUG)
				cout << "Bottom ----> ratio for CH#" << 2*i << "-CH#" << 2*i+1 << " " << ratioBottom[i] << " " << weightedBottom[i] << " " << endl;
			//
			if (data->ADC(2*i+16) == 0 || data->ADC(2*i+1+16) ==0)
			{
				ratioTop[i]    = 0;
				weightedTop[i] = 0;
			}else{
				ratioTop[i] = (Double_t)data->ADC(2*i+16)/data->ADC(2*i+1+16);
				weightedTop[i] = ratioTop[i]*data->ADC(2*i+1+16);
				if(ratioTop[i]>1.) {
					ratioTop[i] = (Double_t)data->ADC(2*i+1+16)/data->ADC(2*i+16);
					weightedTop[i] = ratioTop[i]*data->ADC(2*i+16);
				}
			}
			if(DEBUG) 
				cout << "Top ----> ratio for CH#" << 2*i+16 << "-CH#" << 2*i+1+16 << " " << ratioTop[i] << " " << weightedTop[i] << " " << endl;
		}
		PossibleCoincidence = TMath::LocMax(8,weightedBottom);
		MaxAdcInPosCoincidence = TMath::Max(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
		MinAdcInPosCoincidence = TMath::Min(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
		if(DEBUG) 
		cout << ientry 
			 << " locate max. " 
			 << PossibleCoincidence 
			 << " " << weightedBottom[PossibleCoincidence] 
			 << " CH#" << 2*PossibleCoincidence 
			 << "-CH#" << 2*PossibleCoincidence+1 
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence+1)
			 << " maxADC " << MaxAdcInPosCoincidence
			 << " minADC " << MinAdcInPosCoincidence
			 << endl;
		
		//ScintillatorY[ientry] = PossibleCoincidence;
		
		Long64_t TwoChannelsInPosCoincidence_Plane1 = (2*PossibleCoincidence+2*PossibleCoincidence+1)/2;
		if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane1 " << TwoChannelsInPosCoincidence_Plane1 << endl;
		//
		PossibleCoincidence2 = TMath::LocMax(8,weightedTop);
		MaxAdcInPosCoincidence2 = TMath::Max(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
		MinAdcInPosCoincidence2 = TMath::Min(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
		if(DEBUG) 
		cout << ientry 
			 << " locate max2. " 
			 << PossibleCoincidence2 
			 << " " << weightedTop[PossibleCoincidence2] 
			 << " CH#" << 2*PossibleCoincidence2+16
			 << "-CH#" << 2*PossibleCoincidence2+1+16
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence2+16)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence2+1+16)
			 << " maxADC " << MaxAdcInPosCoincidence2
			 << " minADC " << MinAdcInPosCoincidence2
			 << endl;

		hCoinPlane->Fill(PossibleCoincidence2+1,PossibleCoincidence+1);
		
		//ScintillatorX[ientry] = PossibleCoincidence2;
		Long64_t TwoChannelsInPosCoincidence_Plane2 = (2*PossibleCoincidence2+16+2*PossibleCoincidence2+1+16)/2;
		if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane2 " << TwoChannelsInPosCoincidence_Plane2 << endl;
			for (int i = 0; i < 32; i++ )
		{
			if(DEBUG) cout << i << " " << data->ADC(i) << " " << data->Coincidence() << " " << endl;
			//if(i==31) InCoincidence = ConvertToBinary(ientry,data->Coincidence(),16);
			if(i==31) v_InCoincidence = ConvertToBinary2nd(ientry,data->Coincidence(),16);
			histos->SetBinContent(i+1,data->ADC(i));
		}
		// Count number of coincidence more than one bar!!!
		if(v_InCoincidence.size()>1) multip++;
		
		for(int i=0; i<v_InCoincidence.size();i++)
		{
			if(DEBUG) cout << "Coincidence " << v_InCoincidence[i] << " " << PossibleCoincidence << " " << PossibleCoincidence2 << endl;
			if(DEBUG) cout << "CAEN coincidence pl " << v_InCoincidence[i] << " " << 2*v_InCoincidence[i] << " " << 2*v_InCoincidence[i]+1 << endl;
			h_CAENCoin->SetBinContent(2*v_InCoincidence[i]-1,data->ADC(2*v_InCoincidence[i]-2));
			h_CAENCoin->SetBinContent(2*v_InCoincidence[i],data->ADC(2*v_InCoincidence[i]-1));
		}

		h_PosCoinBottom->SetBinContent(TwoChannelsInPosCoincidence_Plane1+1,data->ADC(TwoChannelsInPosCoincidence_Plane1));
		h_PosCoinBottom->SetBinContent(TwoChannelsInPosCoincidence_Plane1+2,data->ADC(TwoChannelsInPosCoincidence_Plane1+1));
		h_PosCoinTop->SetBinContent(TwoChannelsInPosCoincidence_Plane2+1,data->ADC(TwoChannelsInPosCoincidence_Plane2));
		h_PosCoinTop->SetBinContent(TwoChannelsInPosCoincidence_Plane2+2,data->ADC(TwoChannelsInPosCoincidence_Plane2+1));
		//
		histos->SetMaximum(5000);
		histos->Draw();
		histos->GetYaxis()->SetTitleOffset(1.3);
		h_CAENCoin->SetFillColor(kRed);
		h_CAENCoin->SetFillStyle(1001);
		h_CAENCoin->Draw("same");
		h_PosCoinBottom->SetFillColor(kGreen);
		h_PosCoinBottom->SetFillStyle(3012);
		//h_PosCoinBottom->Draw("same");
		h_PosCoinTop->SetFillColor(kGreen);
		h_PosCoinTop->SetFillStyle(3022);
		//h_PosCoinTop->Draw("same");
		TLine *line = new TLine(15.5,0,15.5,5000);
		line->SetLineColor(kBlue);
		line->SetLineWidth(2);
		line->Draw("same");
		TText *t = new TText(7.5,4500,"Bottom Layer (15mm Bars)");
		t->SetTextAlign(22);
		t->SetTextColor(kBlue);
		t->SetTextFont(43);
		t->SetTextSize(14);
		t->Draw("same");
		TText *t2 = new TText(23,4500,"Top Layer (10mm Bars)");
		t2->SetTextAlign(22);
		t2->SetTextColor(kBlue);
		t2->SetTextFont(43);
		t2->SetTextSize(14);
		t2->Draw("same");
		myCan->Update();
		ps->NewPage();
		//myCan->Clear();
		TH1 *frame = myCan->DrawFrame(0,0,1850,1850);
		frame->GetXaxis()->SetBinLabel(80,"Scin8");
		frame->GetXaxis()->SetBinLabel(200,"Scin7");
		frame->GetXaxis()->SetBinLabel(320,"Scin6");
		frame->GetXaxis()->SetBinLabel(440,"Scin5");
		frame->GetXaxis()->SetBinLabel(560,"Scin4");
		frame->GetXaxis()->SetBinLabel(680,"Scin3");
		frame->GetXaxis()->SetBinLabel(800,"Scin2");
		frame->GetXaxis()->SetBinLabel(920,"Scin1");
		frame->GetXaxis()->SetLabelSize(0.04);

		frame->GetYaxis()->Set(8,frame->GetMinimum(),frame->GetMaximum());
		frame->GetYaxis()->SetBinLabel(8,"Scin16");
		frame->GetYaxis()->SetBinLabel(7,"Scin15");
		frame->GetYaxis()->SetBinLabel(6,"Scin14");
		frame->GetYaxis()->SetBinLabel(5,"Scin13");
		frame->GetYaxis()->SetBinLabel(4,"Scin12");
		frame->GetYaxis()->SetBinLabel(3,"Scin11");
		frame->GetYaxis()->SetBinLabel(2,"Scin10");
		frame->GetYaxis()->SetBinLabel(1,"Scin9");
		frame->GetYaxis()->SetLabelSize(0.04);

		frame->GetXaxis()->SetTitleOffset(1.6);
		frame->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
		frame->GetYaxis()->SetTitleOffset(1.6);
		frame->GetYaxis()->SetTitle("Top Layer (10 mm)");
		frame->SetTitle("CRT Module (Scintillators in Coincidence)");
		myCan->Update();
		//
		LayerY = new TBox(1850-230*PossibleCoincidence,0,1850-230*(PossibleCoincidence+1),1850);
		LayerY->SetFillColor(kGray+1);
		LayerY->SetFillStyle(3006);
		LayerY->SetLineColor(kBlue);
		LayerY->SetLineWidth(2);
		LayerY->Draw();
		LayerX = new TBox(0,230*PossibleCoincidence2,1850,230*PossibleCoincidence2+230);
		LayerX->SetFillColor(kGray+1);
		LayerX->SetFillStyle(3007);
		LayerX->SetLineColor(kRed);
		LayerX->SetLineWidth(2);
		LayerX->Draw("same");
		myCan->Update();
		ps->NewPage();
	}
	hCoinPlane->Draw();
	myCan->Update();
	ps->NewPage();
	ps->Close();
	
	cout << "##############################################" << endl;
	cout << "##############################################" << endl;
	cout << "# Total #of cosmic muon event                #" << endl;
	cout << "    " << total << endl;
	cout << "# Number of bars == 1 Triggering coincidence #" << endl;
	cout << "    " << total-multip << endl;
	cout << "# Number of bars >2 Triggering coincidence   #" << endl;
	cout << "    " << multip << endl;
	cout << "##############################################" << endl;
	cout << "##############################################" << endl;
	
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// get wrongly flagged events for fixed scintillator bars in X-Y
void CheckCoincidenceAllWithoutPlotting(){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	cout << "Number of entries " << entries << endl;
	UInt_t BarX, BarY;
	Int_t num;
	Int_t InCoincidence;
	std::vector<int> v_InCoincidence;
	UInt_t PossibleCoincidence;
	UInt_t PossibleCoincidence2;
	UInt_t MaxAdcInPosCoincidence;
	UInt_t MinAdcInPosCoincidence;
	UInt_t MaxAdcInPosCoincidence2;
	UInt_t MinAdcInPosCoincidence2;
	
	Int_t total = 0, multip = 0, zero = 0, one = 0, two = 0;
	Int_t top_flag=0, bottom_flag=0;
	Double_t ratioBottom[8];
	Double_t weightedBottom[8];
	Double_t ratioTop[8];
	Double_t weightedTop[8];
	Double_t ratio_caen;
	Double_t weighted_caen;
	
	int wrongly = 0, wrongly_b =0, wrongly_t =0;
	TBox *LayerX, *LayerY;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		//if(ientry>100) break;
		data->GetMppc(ientry);
		total++;
		for (int i=0;i<8;i++)
		{
			//cout << "even " << 2*i << " " << event1_adc[2*i] << endl;
			//cout << "odd " << 2*i+1 << " "<< event1_adc[2*i+1] << endl;
			if (data->ADC(2*i) == 0 || data->ADC(2*i+1) ==0)
			{
				ratioBottom[i]    = 0;
				weightedBottom[i] = 0;
			}else{
				ratioBottom[i] = (Double_t)data->ADC(2*i)/data->ADC(2*i+1);
				if(ratioBottom[i]>1.) {
					ratioBottom[i] = (Double_t)data->ADC(2*i+1)/data->ADC(2*i);
				}
				weightedBottom[i] = (ratioBottom[i]+1)*(data->ADC(2*i+1)+data->ADC(2*i));
			}
			if(DEBUG)
				cout << "Bottom ----> ratio for CH#" << 2*i << "-CH#" << 2*i+1 << " " << ratioBottom[i] << " " << weightedBottom[i] << " " << endl;
			//
			if (data->ADC(2*i+16) == 0 || data->ADC(2*i+1+16) ==0)
			{
				ratioTop[i]    = 0;
				weightedTop[i] = 0;
			}else{
				ratioTop[i] = (Double_t)data->ADC(2*i+16)/data->ADC(2*i+1+16);
				if(ratioTop[i]>1.) {
					ratioTop[i] = (Double_t)data->ADC(2*i+1+16)/data->ADC(2*i+16);
				}
				weightedTop[i] = (ratioTop[i]+1)*(data->ADC(2*i+1+16)+data->ADC(2*i+16));
			}
			if(DEBUG) 
				cout << "Top ----> ratio for CH#" << 2*i+16 << "-CH#" << 2*i+1+16 << " " << ratioTop[i] << " " << weightedTop[i] << " " << endl;
		}
		PossibleCoincidence = TMath::LocMax(8,weightedBottom);
		MaxAdcInPosCoincidence = TMath::Max(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
		MinAdcInPosCoincidence = TMath::Min(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
		if(DEBUG) 
		cout << ientry 
			 << " locate max. " 
			 << PossibleCoincidence 
			 << " " << weightedBottom[PossibleCoincidence] 
			 << " CH#" << 2*PossibleCoincidence 
			 << "-CH#" << 2*PossibleCoincidence+1 
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence+1)
			 << " maxADC " << MaxAdcInPosCoincidence
			 << " minADC " << MinAdcInPosCoincidence
			 << endl;
		
		//ScintillatorY[ientry] = PossibleCoincidence;
		
		Long64_t TwoChannelsInPosCoincidence_Plane1 = (2*PossibleCoincidence+2*PossibleCoincidence+1)/2;
		if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane1 " << TwoChannelsInPosCoincidence_Plane1 << endl;
		//
		PossibleCoincidence2 = TMath::LocMax(8,weightedTop);
		MaxAdcInPosCoincidence2 = TMath::Max(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
		MinAdcInPosCoincidence2 = TMath::Min(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
		if(DEBUG) 
		cout << ientry 
			 << " locate max2. " 
			 << PossibleCoincidence2 
			 << " " << weightedTop[PossibleCoincidence2] 
			 << " CH#" << 2*PossibleCoincidence2+16
			 << "-CH#" << 2*PossibleCoincidence2+1+16
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence2+16)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence2+1+16)
			 << " maxADC " << MaxAdcInPosCoincidence2
			 << " minADC " << MinAdcInPosCoincidence2
			 << endl;

		//ScintillatorX[ientry] = PossibleCoincidence2;
		Long64_t TwoChannelsInPosCoincidence_Plane2 = (2*PossibleCoincidence2+16+2*PossibleCoincidence2+1+16)/2;
		if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane2 " << TwoChannelsInPosCoincidence_Plane2 << endl;
			for (int i = 0; i < 32; i++ )
		{
			if(DEBUG) cout << i << " " << data->ADC(i) << " " << data->Coincidence() << " " << endl;
			//if(i==31) InCoincidence = ConvertToBinary(ientry,data->Coincidence(),16);
			if(i==31) v_InCoincidence = ConvertToBinary2nd(ientry,data->Coincidence(),16);
		}
		// Count number of coincidence more than one bar!!!
		if(v_InCoincidence.size()>1) multip++;
		if(v_InCoincidence.size()==0) zero++;
		if(v_InCoincidence.size()==1) one++;
		if(v_InCoincidence.size()==2) two++;
		
		for(int i=0; i<v_InCoincidence.size();i++)
		{
			if(DEBUG) cout << "Coincidence " << v_InCoincidence[i] << " " << PossibleCoincidence << " " << PossibleCoincidence2 << endl;
			if(DEBUG) cout << "CAEN coincidence pl " << v_InCoincidence[i] << " " << 2*v_InCoincidence[i] << " " << 2*v_InCoincidence[i]+1 << endl;
			if((v_InCoincidence[i]<9&&v_InCoincidence[i]!=PossibleCoincidence+1)||(v_InCoincidence[i]>9&&v_InCoincidence[i]-8!=PossibleCoincidence2+1))
				if(DEBUG) cout << "suspicious " << ientry << " " << v_InCoincidence[i] << " " << v_InCoincidence[i]-8 << " " << PossibleCoincidence+1 << " " << PossibleCoincidence2+1 << endl;
		}
	
		if(v_InCoincidence.size()==1)
		{
			cout << ientry << " CAEN coincidence pl " << v_InCoincidence[0] << " " << 2*v_InCoincidence[0]-1 << " " << 2*v_InCoincidence[0]-2 << endl;
			cout << ientry << " " << data->ADC(2*v_InCoincidence[0]-1) << " " << data->ADC(2*v_InCoincidence[0]-2) << endl;			
			
			if(data->ADC(2*v_InCoincidence[0]-1)==0||data->ADC(2*v_InCoincidence[0]-2)==0)
			{
				cout << "umut" << endl;
				ratio_caen = 0;
			} else {
				ratio_caen = (Double_t) data->ADC(2*v_InCoincidence[0]-1)/data->ADC(2*v_InCoincidence[0]-2);
				if(ratio_caen>1.)
				{
					ratio_caen = (Double_t) data->ADC(2*v_InCoincidence[0]-2)/data->ADC(2*v_InCoincidence[0]-1);
					cout << "greater ratio" << endl;
				}
				weighted_caen = (ratio_caen+1)*(data->ADC(2*v_InCoincidence[0]-1) + data->ADC(2*v_InCoincidence[0]-2));
				cout << " .. " << ratio_caen << " " << weighted_caen << " " << data->ADC(2*v_InCoincidence[0]-1) << " " << data->ADC(2*v_InCoincidence[0]-2) << endl;
			}
			if(v_InCoincidence[0]>8)
			{
				cout << ientry << " Top layer " << endl;
				if(weightedBottom[PossibleCoincidence] < weighted_caen) 
				{
					cout << ientry << " wrongly top flaged " << endl;
					wrongly_t++;
					wrongly++;
				}
			}
			if(v_InCoincidence[0]<9)
			{
				cout << ientry << " Bottom layer " << endl;
				if(weightedTop[PossibleCoincidence2] < weighted_caen) 
				{
					cout << ientry << " wrongly bottom flaged " << endl;
					wrongly_b++;
					wrongly++;
				}
			}
			cout << ientry << " caen: " << weighted_caen << " bot: " << weightedBottom[PossibleCoincidence] << " top: " << weightedTop[PossibleCoincidence2] << endl;
		}
		if(weightedTop[PossibleCoincidence2] < weightedBottom[PossibleCoincidence]) 
		{
			cout << ientry << " Top layer should be flaged " << endl;
			top_flag++;
		}
		if(weightedTop[PossibleCoincidence2] > weightedBottom[PossibleCoincidence]) 
		{
			cout << ientry << " Bottom layer should be flaged " << endl;
			bottom_flag++;
		}
	}
	
	
	cout << "##############################################" << endl;
	cout << "##############################################" << endl;
	cout << "# Total #of cosmic muon event                #" << endl;
	cout << "    " << total << endl;
	cout << "# Number of bars == 0 Triggering coincidence #" << endl;
	cout << "    " << zero << " " << zero/total << endl;
	cout << "# Number of bars == 1 Triggering coincidence #" << endl;
	cout << "    " << one << " " << one/total << endl;
	cout << 	total << " " << wrongly_t << " " << wrongly_b << " " << wrongly << endl;
	cout << "# Number of bars == 2 Triggering coincidence   #" << endl;
	cout << "    " << two << " " << two/total << endl;
	cout << "# Number of bars > 1 Triggering coincidence   #" << endl;
	cout << "    " << multip << " " << multip/total << endl;
	cout << "##############################################" << endl;
	cout << "##############################################" << endl;
	cout << "Flag should be top    " << top_flag << endl;
	cout << "Flag should be bottom " << bottom_flag << endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void CheckCoincidenceSingle(int ientry){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
//	TPostScript *ps = new TPostScript(Form("Coincidence_entry%i.ps",ientry),112);
	TH1F *m = (TH1F*) gROOT->FindObject("histos");
	if(m) delete m;
	m = (TH1F*) gROOT->FindObject("h_CAENCoin");
	if(m) delete m;
	m = (TH1F*) gROOT->FindObject("h_PosCoinBottom");
	if(m) delete m;
	m = (TH1F*) gROOT->FindObject("h_PosCoinTop");
	if(m) delete m;
	TH2F *m2 = (TH2F*) gROOT->FindObject("hCoinPlane");
	if(m2) delete m2;

	TH1F *histos = new TH1F("histos",Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);;
	TH1F *h_CAENCoin = new TH1F("h_CAENCoin",Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);;
	TH1F *h_PosCoinBottom = new TH1F("h_PosCoinBottom",Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);;
	TH1F *h_PosCoinTop = new TH1F("h_PosCoinTop",Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);;
	Int_t ScintillatorX;
	Int_t ScintillatorY;
	
	TH2F *hCoinPlane = new TH2F("hCoinPlane","Cosmic Muon Event in Coincidence;15 mm Scintillator Bar;10 mm Scintillator Bar",20,0,8,20,0,8);

//	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
//	if (C) delete C;
//	TCanvas *myCan = new TCanvas("Canvas","Quality Check",200,10,700,500);

	Mppc *data = filereader.tmCD();
	UInt_t BarX, BarY;
	Int_t num;
	Int_t InCoincidence;
	std::vector<int> v_InCoincidence;
	UInt_t PossibleCoincidence;
	UInt_t PossibleCoincidence2;
	UInt_t MaxAdcInPosCoincidence;
	UInt_t MinAdcInPosCoincidence;
	UInt_t MaxAdcInPosCoincidence2;
	UInt_t MinAdcInPosCoincidence2;

	Double_t ratioBottom[8];
	Double_t weightedBottom[8];
	Double_t ratioTop[8];
	Double_t weightedTop[8];
	TBox *LayerX, *LayerY;
	//
	data->GetMppc(ientry);
	for (int i=0;i<8;i++)
	{
		if (data->ADC(2*i) == 0 || data->ADC(2*i+1) ==0)
		{
			ratioBottom[i]    = 0;
			weightedBottom[i] = 0;
		}else{
			ratioBottom[i] = (Double_t)data->ADC(2*i)/data->ADC(2*i+1);
			weightedBottom[i] = (ratioBottom[i]+1)*(data->ADC(2*i+1)+data->ADC(2*i));
			if(ratioBottom[i]>1.) {
				ratioBottom[i] = (Double_t)data->ADC(2*i+1)/data->ADC(2*i);
				weightedBottom[i] = (ratioBottom[i]+1)*(data->ADC(2*i)+data->ADC(2*i+1));
			}
		}
		if(DEBUG1) 
			cout << "Bottom ----> ratio for CH#" << 2*i << "-CH#" << 2*i+1 << " " << ratioBottom[i] << " " << weightedBottom[i] << " " << endl;
		//
		if (data->ADC(2*i+16) == 0 || data->ADC(2*i+1+16) ==0)
		{
			ratioTop[i]    = 0;
			weightedTop[i] = 0;
		}else{
			ratioTop[i] = (Double_t)data->ADC(2*i+16)/data->ADC(2*i+1+16);
			weightedTop[i] = (ratioTop[i]+1)*(data->ADC(2*i+1+16)+data->ADC(2*i+16));
		if(ratioTop[i]>1.) {
				ratioTop[i] = (Double_t)data->ADC(2*i+1+16)/data->ADC(2*i+16);
				weightedTop[i] = (ratioTop[i]+1)*(data->ADC(2*i+16)+data->ADC(2*i+1+16));
			}
		}
		if(DEBUG1) 
			cout << "Top ----> ratio for CH#" << 2*i+16 << "-CH#" << 2*i+1+16 << " " << ratioTop[i] << " " << weightedTop[i] << " " << endl;
	}
	PossibleCoincidence = TMath::LocMax(8,weightedBottom);
	MaxAdcInPosCoincidence = TMath::Max(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
	MinAdcInPosCoincidence = TMath::Min(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
	if(DEBUG1) 
	cout << ientry 
		 << " locate max. " 
		 << PossibleCoincidence 
		 << " " << weightedBottom[PossibleCoincidence] 
		 << " CH#" << 2*PossibleCoincidence 
		 << "-CH#" << 2*PossibleCoincidence+1 
		 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence)
		 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence+1)
		 << " maxADC " << MaxAdcInPosCoincidence
		 << " minADC " << MinAdcInPosCoincidence
		 << endl;
	
	ScintillatorY = PossibleCoincidence;
	
	Long64_t TwoChannelsInPosCoincidence_Plane1 = (2*PossibleCoincidence+2*PossibleCoincidence+1)/2;
	if(DEBUG1) cout << "TwoChannelsInPosCoincidence_Plane1 " << TwoChannelsInPosCoincidence_Plane1 << endl;
	//
	PossibleCoincidence2 = TMath::LocMax(8,weightedTop);
	MaxAdcInPosCoincidence2 = TMath::Max(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
	MinAdcInPosCoincidence2 = TMath::Min(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
	if(DEBUG1) 
	cout << ientry 
		 << " locate max2. " 
		 << PossibleCoincidence2 
		 << " " << weightedTop[PossibleCoincidence2] 
		 << " CH#" << 2*PossibleCoincidence2+16
		 << "-CH#" << 2*PossibleCoincidence2+1+16
		 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence2+16)
		 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence2+1+16)
		 << " maxADC " << MaxAdcInPosCoincidence2
		 << " minADC " << MinAdcInPosCoincidence2
		 << endl;

	hCoinPlane->Fill(PossibleCoincidence2+1,PossibleCoincidence+1);
		
	ScintillatorX = PossibleCoincidence2;
	Long64_t TwoChannelsInPosCoincidence_Plane2 = (2*PossibleCoincidence2+16+2*PossibleCoincidence2+1+16)/2;
	if(DEBUG1) cout << "TwoChannelsInPosCoincidence_Plane2 " << TwoChannelsInPosCoincidence_Plane2 << endl;
	for (int i = 0; i < 32; i++ )
	{
		if(DEBUG1) cout << i << " " << data->ADC(i) << " " << data->Coincidence() << " " << endl;
		//if(i==31) InCoincidence = ConvertToBinary(ientry,data->Coincidence(),16);
		if(i==31) v_InCoincidence = ConvertToBinary2nd(ientry,data->Coincidence(),16);
		if(DEBUG1) cout << endl;
		histos->SetBinContent(i+1,data->ADC(i));
	}

	//if(DEBUG) cout << "Coincidence " << InCoincidence << " " << PossibleCoincidence << " " << PossibleCoincidence2 << endl;
	//if(DEBUG) cout << "CAEN coincidence pl " << InCoincidence << " " << 2*InCoincidence << " " << 2*InCoincidence+1 << endl;
	//h_CAENCoin->SetBinContent(2*InCoincidence-1,data->ADC(2*InCoincidence-2));
	//h_CAENCoin->SetBinContent(2*InCoincidence,data->ADC(2*InCoincidence-1));

	for(int i=0; i<v_InCoincidence.size();i++)
	{
		if(DEBUG1) cout << "Coincidence " << v_InCoincidence[i] << " " << PossibleCoincidence << " " << PossibleCoincidence2 << endl;
		if(DEBUG1) cout << "CAEN coincidence pl " << v_InCoincidence[i] << " " << 2*v_InCoincidence[i] << " " << 2*v_InCoincidence[i]+1 << endl;
		h_CAENCoin->SetBinContent(2*v_InCoincidence[i]-1,data->ADC(2*v_InCoincidence[i]-2));
		h_CAENCoin->SetBinContent(2*v_InCoincidence[i],data->ADC(2*v_InCoincidence[i]-1));
	}
	
	
	h_PosCoinBottom->SetBinContent(TwoChannelsInPosCoincidence_Plane1+1,data->ADC(TwoChannelsInPosCoincidence_Plane1));
	h_PosCoinBottom->SetBinContent(TwoChannelsInPosCoincidence_Plane1+2,data->ADC(TwoChannelsInPosCoincidence_Plane1+1));
	h_PosCoinTop->SetBinContent(TwoChannelsInPosCoincidence_Plane2+1,data->ADC(TwoChannelsInPosCoincidence_Plane2));
	h_PosCoinTop->SetBinContent(TwoChannelsInPosCoincidence_Plane2+2,data->ADC(TwoChannelsInPosCoincidence_Plane2+1));

	c1->cd();
	histos->SetMaximum(5000);
	histos->Draw();
	histos->GetYaxis()->SetTitleOffset(1.3);
	h_CAENCoin->SetFillColor(kRed);
	h_CAENCoin->SetFillStyle(1001);
	h_CAENCoin->Draw("same");
	h_PosCoinBottom->SetFillColor(kGreen);
	h_PosCoinBottom->SetFillStyle(3012);
	h_PosCoinBottom->Draw("same");
	h_PosCoinTop->SetFillColor(kGreen);
	h_PosCoinTop->SetFillStyle(3022);
	h_PosCoinTop->Draw("same");
	TLine *line = new TLine(15.5,0,15.5,5000);
	line->SetLineColor(kBlue);
	line->SetLineWidth(2);
	line->Draw("same");
	TText *t = new TText(7.5,4500,"Bottom Layer (15mm Bars)");
	t->SetTextAlign(22);
	t->SetTextColor(kBlue);
	t->SetTextFont(43);
	t->SetTextSize(14);
	t->Draw("same");
	TText *t2 = new TText(23,4500,"Top Layer (10mm Bars)");
	t2->SetTextAlign(22);
	t2->SetTextColor(kBlue);
	t2->SetTextFont(43);
	t2->SetTextSize(14);
	t2->Draw("same");
	c1->Modified();
	c1->Update();
//	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PlotInCoincidence(){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TPostScript *ps = new TPostScript("InCoincidence.ps",112);
	
	TH2F *hCoinPlane = new TH2F("hCoinPlane","Cosmic Muon Event in Coincidence;15 mm Scintillator Bar;10 mm Scintillator Bar",20,0,9,20,0,9);
	TH2F *hCoinPlane2nd = new TH2F("hCoinPlane2nd","Cosmic Muon Event in Coincidence;15 mm Scintillator Bar;10 mm Scintillator Bar",20,0,9,20,0,9);
	//
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","In Coincidence",200,10,700,500);
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	cout << "Number of entries " << entries << endl;
	UInt_t BarX, BarY;
	Int_t num;
	Int_t InCoincidence;
	std::vector<int> v_InCoincidence;
	UInt_t PossibleCoincidence;
	UInt_t PossibleCoincidence2;
	UInt_t MaxAdcInPosCoincidence;
	UInt_t MinAdcInPosCoincidence;
	UInt_t MaxAdcInPosCoincidence2;
	UInt_t MinAdcInPosCoincidence2;

	Double_t ratioBottom[8];
	Double_t weightedBottom[8];
	Double_t ratioTop[8];
	Double_t weightedTop[8];
	TBox *LayerX, *LayerY;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		//if(ientry >10000) continue;
		data->GetMppc(ientry);
		for (int i=0;i<8;i++)
		{
			//cout << "even " << 2*i << " " << event1_adc[2*i] << endl;
			//cout << "odd " << 2*i+1 << " "<< event1_adc[2*i+1] << endl;
			if (data->ADC(2*i) == 0 || data->ADC(2*i+1) ==0)
			{
				ratioBottom[i]    = 0;
				weightedBottom[i] = 0;
			}else{
				ratioBottom[i] = (Double_t)data->ADC(2*i)/data->ADC(2*i+1);
				weightedBottom[i] = ratioBottom[i]*data->ADC(2*i+1);
				if(ratioBottom[i]>1.) {
					ratioBottom[i] = (Double_t)data->ADC(2*i+1)/data->ADC(2*i);
					weightedBottom[i] = ratioBottom[i]*data->ADC(2*i);
				}
			}
			if(DEBUG) 
				cout << "Bottom ----> ratio for CH#" << 2*i << "-CH#" << 2*i+1 << " " << ratioBottom[i] << " " << weightedBottom[i] << " " << endl;
			//
			if (data->ADC(2*i+16) == 0 || data->ADC(2*i+1+16) ==0)
			{
				ratioTop[i]    = 0;
				weightedTop[i] = 0;
			}else{
				ratioTop[i] = (Double_t)data->ADC(2*i+16)/data->ADC(2*i+1+16);
				weightedTop[i] = ratioTop[i]*data->ADC(2*i+1+16);
				if(ratioTop[i]>1.) {
					ratioTop[i] = (Double_t)data->ADC(2*i+1+16)/data->ADC(2*i+16);
					weightedTop[i] = ratioTop[i]*data->ADC(2*i+16);
				}
			}
			if(DEBUG) 
				cout << "Top ----> ratio for CH#" << 2*i+16 << "-CH#" << 2*i+1+16 << " " << ratioTop[i] << " " << weightedTop[i] << " " << endl;
		}
		PossibleCoincidence = TMath::LocMax(8,weightedBottom);
		MaxAdcInPosCoincidence = TMath::Max(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
		MinAdcInPosCoincidence = TMath::Min(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
		if(DEBUG) 
		cout << ientry 
			 << " locate max. " 
			 << PossibleCoincidence 
			 << " " << weightedBottom[PossibleCoincidence] 
			 << " CH#" << 2*PossibleCoincidence 
			 << "-CH#" << 2*PossibleCoincidence+1 
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence+1)
			 << " maxADC " << MaxAdcInPosCoincidence
			 << " minADC " << MinAdcInPosCoincidence
			 << endl;
		
		Long64_t TwoChannelsInPosCoincidence_Plane1 = (2*PossibleCoincidence+2*PossibleCoincidence+1)/2;
		if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane1 " << TwoChannelsInPosCoincidence_Plane1 << endl;
		//
		PossibleCoincidence2 = TMath::LocMax(8,weightedTop);
		MaxAdcInPosCoincidence2 = TMath::Max(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
		MinAdcInPosCoincidence2 = TMath::Min(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
		if(DEBUG) 
		cout << ientry 
			 << " locate max2. " 
			 << PossibleCoincidence2 
			 << " " << weightedTop[PossibleCoincidence2] 
			 << " CH#" << 2*PossibleCoincidence2+16
			 << "-CH#" << 2*PossibleCoincidence2+1+16
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence2+16)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence2+1+16)
			 << " maxADC " << MaxAdcInPosCoincidence2
			 << " minADC " << MinAdcInPosCoincidence2
			 << endl;

		for (int i = 0; i < 32; i++ )
		{
			if(i==31) InCoincidence = ConvertToBinary(ientry,data->Coincidence(),16);
			if(i==31) v_InCoincidence = ConvertToBinary2nd(ientry,data->Coincidence(),16);
		}
		hCoinPlane->Fill(PossibleCoincidence+1,PossibleCoincidence2+1);
		//
		for(int i=0; i<v_InCoincidence.size();i++)
		{
			if(v_InCoincidence[i]!=PossibleCoincidence2+9) // top layer bars: 9-10-11-12-13-14-15-16
			hCoinPlane2nd->Fill(v_InCoincidence[i]-8,PossibleCoincidence2+1);
			if(v_InCoincidence[i]!=PossibleCoincidence+1) // bottom layer bars: 1-2-3-4-5-6-7-8
			hCoinPlane2nd->Fill(PossibleCoincidence+1,v_InCoincidence[i]);
		}
		//
		Long64_t TwoChannelsInPosCoincidence_Plane2 = (2*PossibleCoincidence2+16+2*PossibleCoincidence2+1+16)/2;
		if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane2 " << TwoChannelsInPosCoincidence_Plane2 << endl;
	}
	myCan->SetGrid();
	hCoinPlane->Draw("text");
	myCan->Update();
	ps->NewPage();
	hCoinPlane2nd->Draw("text");
	myCan->Update();
	ps->NewPage();
	ps->Close();
}








//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PlotInCoincidenceSecond(){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TPostScript *ps = new TPostScript("InCoincidence.ps",112);
	
	TH2F *hCoinPlane = new TH2F("hCoinPlane","Cosmic Muon Event in Coincidence;15 mm Scintillator Bar;10 mm Scintillator Bar",20,0,9,20,0,9);
	//
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","In Coincidence",200,10,700,500);
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	cout << "Number of entries " << entries << endl;
	UInt_t BarX, BarY;
	Int_t num;
	Int_t InCoincidence;
	UInt_t PossibleCoincidence;
	UInt_t PossibleCoincidence2;
	UInt_t MaxAdcInPosCoincidence;
	UInt_t MinAdcInPosCoincidence;
	UInt_t MaxAdcInPosCoincidence2;
	UInt_t MinAdcInPosCoincidence2;

	Double_t ratioBottom[8];
	Double_t weightedBottom[8];
	Double_t ratioTop[8];
	Double_t weightedTop[8];
	TBox *LayerX, *LayerY;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		//if(ientry >10000) continue;
		data->GetMppc(ientry);
		for (int i=0;i<8;i++)
		{
			//cout << "even " << 2*i << " " << event1_adc[2*i] << endl;
			//cout << "odd " << 2*i+1 << " "<< event1_adc[2*i+1] << endl;
			if (data->ADC(2*i) == 0 || data->ADC(2*i+1) ==0)
			{
				ratioBottom[i]    = 0;
				weightedBottom[i] = 0;
			}else{
				ratioBottom[i] = (Double_t)data->ADC(2*i)/data->ADC(2*i+1);
				weightedBottom[i] = ratioBottom[i]*data->ADC(2*i+1);
				if(ratioBottom[i]>1.) {
					ratioBottom[i] = (Double_t)data->ADC(2*i+1)/data->ADC(2*i);
					weightedBottom[i] = ratioBottom[i]*data->ADC(2*i);
				}
			}
			if(DEBUG) 
				cout << "Bottom ----> ratio for CH#" << 2*i << "-CH#" << 2*i+1 << " " << ratioBottom[i] << " " << weightedBottom[i] << " " << endl;
			//
			if (data->ADC(2*i+16) == 0 || data->ADC(2*i+1+16) ==0)
			{
				ratioTop[i]    = 0;
				weightedTop[i] = 0;
			}else{
				ratioTop[i] = (Double_t)data->ADC(2*i+16)/data->ADC(2*i+1+16);
				weightedTop[i] = ratioTop[i]*data->ADC(2*i+1+16);
				if(ratioTop[i]>1.) {
					ratioTop[i] = (Double_t)data->ADC(2*i+1+16)/data->ADC(2*i+16);
					weightedTop[i] = ratioTop[i]*data->ADC(2*i+16);
				}
			}
			if(DEBUG) 
				cout << "Top ----> ratio for CH#" << 2*i+16 << "-CH#" << 2*i+1+16 << " " << ratioTop[i] << " " << weightedTop[i] << " " << endl;
		}
		PossibleCoincidence = TMath::LocMax(8,weightedBottom);
		MaxAdcInPosCoincidence = TMath::Max(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
		MinAdcInPosCoincidence = TMath::Min(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
		if(DEBUG) 
		cout << ientry 
			 << " locate max. " 
			 << PossibleCoincidence 
			 << " " << weightedBottom[PossibleCoincidence] 
			 << " CH#" << 2*PossibleCoincidence 
			 << "-CH#" << 2*PossibleCoincidence+1 
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence+1)
			 << " maxADC " << MaxAdcInPosCoincidence
			 << " minADC " << MinAdcInPosCoincidence
			 << endl;
		
		Long64_t TwoChannelsInPosCoincidence_Plane1 = (2*PossibleCoincidence+2*PossibleCoincidence+1)/2;
		if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane1 " << TwoChannelsInPosCoincidence_Plane1 << endl;
		//
		PossibleCoincidence2 = TMath::LocMax(8,weightedTop);
		MaxAdcInPosCoincidence2 = TMath::Max(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
		MinAdcInPosCoincidence2 = TMath::Min(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
		if(DEBUG) 
		cout << ientry 
			 << " locate max2. " 
			 << PossibleCoincidence2 
			 << " " << weightedTop[PossibleCoincidence2] 
			 << " CH#" << 2*PossibleCoincidence2+16
			 << "-CH#" << 2*PossibleCoincidence2+1+16
			 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence2+16)
			 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence2+1+16)
			 << " maxADC " << MaxAdcInPosCoincidence2
			 << " minADC " << MinAdcInPosCoincidence2
			 << endl;

		hCoinPlane->Fill(PossibleCoincidence2+1,PossibleCoincidence+1);
		
		Long64_t TwoChannelsInPosCoincidence_Plane2 = (2*PossibleCoincidence2+16+2*PossibleCoincidence2+1+16)/2;
		if(DEBUG) cout << "TwoChannelsInPosCoincidence_Plane2 " << TwoChannelsInPosCoincidence_Plane2 << endl;
	}
	myCan->SetGrid();
	hCoinPlane->Draw("text");
	myCan->Update();
	ps->NewPage();
	ps->Close();
}








//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#define DEBUG2 1

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int ConvertToBinary(int entry, int number, int num_digits) {
	int digit;
	int q=0;
	int BarNumber;
	int NumOf1s = 0;
	for(digit = num_digits - 1; digit >= 0; digit--)
	{
		if(DEBUG2) printf("%c", number & (1 << digit) ? '1' : '0');
		q++;
		if (q==4)
		{
			if(DEBUG2) printf(" ");
			q=0;
		}
		if(number & (1 << digit)) NumOf1s++;
		if(number & (1 << digit)) BarNumber = digit+1;
	}
	if(DEBUG2) printf("\n");
	if(NumOf1s>1) cout << "Number of 1s more than one for entry " << entry << " !!!!!" << endl;
	if(DEBUG2) printf("%d %d",NumOf1s,BarNumber);
	return BarNumber;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std::vector<int> ConvertToBinary2nd(int entry, int number, int num_digits) {
	int digit;
	int q=0;
	std::vector<int> BarNumber;
	int NumOf1s = 0;
	for(digit = num_digits - 1; digit >= 0; digit--)
	{
		if(DEBUG1) printf("%c", number & (1 << digit) ? '1' : '0');
		q++;
		if (q==4)
		{
			if(DEBUG1) printf(" ");
			q=0;
		}
		if(number & (1 << digit)) NumOf1s++;
		if(number & (1 << digit)) 
		{
			BarNumber.push_back(digit+1);
		}
	}
	if(DEBUG1) printf("\n");
	if(DEBUG1) cout << "Caen coincidence " << BarNumber.size() << " " << entry << endl;
	if(DEBUG1&&NumOf1s>1) cout << "Number of 1s more than one for entry " << entry << " is " << BarNumber.size() << " !!!!!" << endl;
	if(DEBUG1) printf("%d %zu\n",NumOf1s,BarNumber.size());
	return BarNumber;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void TwoFibersInBar(){
	TH2F *hist2[16];
	TH2F *hist2c[16];
	//TPostScript *ps = new TPostScript("TwoFibersInBar.ps",112);
	//ps->NewPage();
	for(int i = 0; i < 8; i++)
	{
		c2->cd(i+1);
		mppc->Draw(Form("chg[%i]:chg[%i] >> htmp%i(1000,0,4100,500,0,4100)",2*i,2*i+1,i));
		hist2[i] = (TH2F*)gDirectory->Get(Form("htmp%i",i));
		hist2c[i] = (TH2F*)hist2[i]->Clone();
		hist2[i]->SetTitle(Form("Two Fibers on 15 mm Scintillator Bar #%i;Right Fiber [CH%i];Left Fiber [CH%i]",(2*i+2)/2,2*i,2*i+1));
		hist2[i]->SetStats(kFALSE);
		hist2[i]->Draw();
		c2->Update();
		if(i==7) c2->Print("TwoFibersInBar.pdf(");
	}
	//ps->NewPage();
	for(int i = 8; i < 16; i++)
	{
		c2b->cd(i-7);
		mppc->Draw(Form("chg[%i]:chg[%i] >> htmp%i(1000,0,4100,500,0,4100)",2*i,2*i+1,i));
		hist2[i] = (TH2F*)gDirectory->Get(Form("htmp%i",i));
		hist2c[i] = (TH2F*)hist2[i]->Clone();
		hist2[i]->SetTitle(Form("Two Fibers on 10 mm Scintillator Bar #%i;Right Fiber [CH%i];Left Fiber [CH%i]",(i-7),2*i,2*i+1));
		hist2[i]->SetStats(kFALSE);
		hist2[i]->Draw();
		c2b->Update();
		if(i==15) c2->Print("TwoFibersInBar.pdf)");
	}
	//ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void TwoFibersInBar_first(){
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Two Fibers in a Bar",200,10,700,500);
	myCan->Divide(2,4);
	myCan->Update();
	myCan->Modified();
	TH2F *hist2[16];
	TH2F *hist2c[16];
	//TPostScript *ps = new TPostScript("TwoFibersInBar.ps",112);
	//ps->NewPage();
	for(int i = 0; i < 16; i++)
	{
		myCan->cd(i+1);
		if(i>7) myCan->cd(i-7);
		mppc->Draw(Form("chg[%i]:chg[%i] >> htmp%i(1000,0,4100,500,0,4100)",2*i,2*i+1,i));
		hist2[i] = (TH2F*)gDirectory->Get(Form("htmp%i",i));
		hist2c[i] = (TH2F*)hist2[i]->Clone();
		if(i < 8 ) hist2[i]->SetTitle(Form("Two Fibers on 15 mm Scintillator Bar #%i;Right Fiber [CH%i];Left Fiber [CH%i]",(2*i+2)/2,2*i,2*i+1));
		if(i > 7 ) hist2[i]->SetTitle(Form("Two Fibers on 10 mm Scintillator Bar #%i;Right Fiber [CH%i];Left Fiber [CH%i]",(i-7),2*i,2*i+1));
		hist2[i]->SetStats(kFALSE);
		hist2[i]->Draw();
		myCan->Update();
		if(i==7) myCan->Print("TwoFibersInBar.pdf(");//ps->NewPage();
		if(i==15) myCan->Print("TwoFibersInBar.pdf)");
	}
	myCan->Update();
	myCan->Clear();
	myCan->Divide(1,1);
	for(int i = 0; i < 16; i++)
	{
		myCan->cd(1);
		if(i < 8 ) hist2c[i]->SetTitle(Form("Two Fibers on 15 mm Scintillator Bar #%i;Right Fiber [CH%i];Left Fiber [CH%i]",(2*i+2)/2,2*i,2*i+1));
		if(i > 7 ) hist2c[i]->SetTitle(Form("Two Fibers on 10 mm Scintillator Bar #%i;Right Fiber [CH%i];Left Fiber [CH%i]",(i-7),2*i,2*i+1));
		hist2c[i]->SetStats(kFALSE);
		//hist2c[i]->SetMarkerStyle(2);
		hist2c[i]->GetYaxis()->SetTitleOffset(1.4);
		hist2c[i]->Draw();
		hist2c[i]->SetMarkerStyle(20);
		hist2c[i]->SetMarkerSize(0.2);
		myCan->Update();
		//ps->NewPage();
	}
	
	
	//ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BarsInTwoPlanes(){
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Bars in Two Planes",200,10,700,500);
	myCan->Clear();
	myCan->Divide(2,4);
	TH2F *hpl2d[8][8];
	TH2F *hpl2dc[8][8];
	int num, index = 0;
//	TPostScript *ps = new TPostScript("BarsInTwoPlanes.ps",112);
//	ps->NewPage();
//	TPDF *ps = new TPDF("BarsInTwoPlanes.pdf");
//	ps->NewPage();
	int j;
	for(int i = 0; i < 8; i++)
	{
		for(j = 8; j < 16; j++)
		{
			num = abs(8-j);
			index++;
			myCan->cd(num+1);
			mppc->Draw(Form("chg[%i]+chg[%i]:chg[%i]+chg[%i] >> htemp%i(500,0,8500,500,0,8500)",2*i,2*i+1,2*j,2*j+1,index));
			hpl2d[i][num] = (TH2F*)gDirectory->Get(Form("htemp%i",index));
			hpl2dc[i][num] = (TH2F*)hpl2d[i][num]->Clone();
			hpl2d[i][num]->SetStats(kFALSE);
			hpl2d[i][num]->SetTitle(Form("15mm Scin#%i vs 10mm Scin#%i;chg[%i]+chg[%i];chg[%i]+chg[%i]",i+1,num+1,2*i,2*i+1,2*j,2*j+1));
			hpl2d[i][num]->Draw();
			myCan->Update();
			myCan->Modified();

		}
		if(i==0) myCan->Print("BarsInTwoPlanes.pdf(");
		if(i!=0) myCan->Print("BarsInTwoPlanes.pdf");
//		ps->NewPage();
	}
	myCan->Update();
	myCan->Clear();
	myCan->Divide(1,1);
	for(int i = 0; i < 8; i++)
	{
		for( j = 8; j < 16; j++)
		{
			myCan->cd(1);
			num = abs(8-j);
			hpl2dc[i][num]->SetStats(kFALSE);
			hpl2dc[i][num]->SetTitle(Form("15mm Scin#%i vs 10mm Scin#%i;chg[%i]+chg[%i];chg[%i]+chg[%i]",i+1,num+1,2*i,2*i+1,2*j,2*j+1));
			hpl2dc[i][num]->Draw();
			hpl2dc[i][num]->SetMarkerStyle(20);
			hpl2dc[i][num]->SetMarkerSize(0.2);
			hpl2dc[i][num]->GetYaxis()->SetTitleOffset(1.4);
			myCan->Update();
			myCan->Modified();
		        myCan->Print("BarsInTwoPlanes.pdf");

//			ps->NewPage();
		}
	}
	myCan->Print("BarsInTwoPlanes.pdf)");
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BarsInTwoPlanes_first(){
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Bars in Two Planes",200,10,700,500);
	myCan->Clear();
	myCan->Divide(2,4);
	TH2F *hpl2d[8][8];
	TH2F *hpl2dc[8][8];
	int num, index = 0;
//	TPostScript *ps = new TPostScript("BarsInTwoPlanes.ps",112);
//	ps->NewPage();
	TPDF *ps = new TPDF("BarsInTwoPlanes.pdf");
//	ps->NewPage();
	for(int i = 0; i < 8; i++)
	{
		for(int j = 8; j < 16; j++)
		{
			num = abs(8-j);
			index++;
			myCan->cd(num+1);
			mppc->Draw(Form("chg[%i]+chg[%i]:chg[%i]+chg[%i] >> htemp%i(500,0,8500,500,0,8500)",2*i,2*i+1,2*j,2*j+1,index));
			hpl2d[i][num] = (TH2F*)gDirectory->Get(Form("htemp%i",index));
			hpl2dc[i][num] = (TH2F*)hpl2d[i][num]->Clone();
			hpl2d[i][num]->SetStats(kFALSE);
			hpl2d[i][num]->SetTitle(Form("15mm Scin#%i vs 10mm Scin#%i;chg[%i]+chg[%i];chg[%i]+chg[%i]",i+1,num+1,2*i,2*i+1,2*j,2*j+1));
			hpl2d[i][num]->Draw();
			myCan->Update();
		}
//		ps->NewPage();
	}
	myCan->Update();
	myCan->Clear();
	myCan->Divide(1,1);
	for(int i = 0; i < 8; i++)
	{
		for(int j = 8; j < 16; j++)
		{
			myCan->cd(1);
			num = abs(8-j);
			hpl2dc[i][num]->SetStats(kFALSE);
			hpl2dc[i][num]->SetTitle(Form("15mm Scin#%i vs 10mm Scin#%i;chg[%i]+chg[%i];chg[%i]+chg[%i]",i+1,num+1,2*i,2*i+1,2*j,2*j+1));
			hpl2dc[i][num]->Draw();
			hpl2dc[i][num]->SetMarkerStyle(20);
			hpl2dc[i][num]->SetMarkerSize(0.2);
			hpl2dc[i][num]->GetYaxis()->SetTitleOffset(1.4);
			myCan->Update();
//			ps->NewPage();
		}
	}
	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void BarsInTwoPlanesDump(){
	TH2F *hpl2d[8][8];
	for(int i = 0; i < 8; i++)
	{
		for(int j = 0; j < 8; j++)
		{
			const char* name = Form("hpl2d_%i_%i",i,j);
			hpl2d[i][j] = new TH2F(name,name,1000,0,8100,1000,0,8100);
		}
	}
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	cout << "Number of entries " << entries << endl;
	UInt_t BarX, BarY;
	int num;
	for( Long64_t ientry=0; ientry<entries; ++ientry )
	{
		data->GetMppc(ientry);
		//
		for(int i = 0; i < 8; i++)
		{
			for(int j = 8; j < 16; j++)
			{
				BarX = data->ADC(2*i) + data->ADC(2*i+1);
				BarY = data->ADC(2*j) + data->ADC(2*j+1);
				//cout << 2*i << " " << 2*i+1 << " " << 2
				num = abs(8-j);
				hpl2d[i][num]->Fill(BarX,BarY);
				cout << ientry << " " << i << " " << num << " " << BarX << " " << BarY << endl;
			}
		}
	}

	TCanvas *myCan = new TCanvas();
	myCan->SaveAs("15mmBarTo10mmBar.pdf(");
	for(int i = 0; i < 8; i++)
	{
		for(int j = 0; j < 8; j++)
		{
			hpl2d[i][j]->SetLineWidth(2);
			hpl2d[i][j]->SetMarkerStyle(2);
			hpl2d[i][j]->SetStats(kFALSE);
			hpl2d[i][j]->SetTitle(Form("15 mm Bar-%i vs 10 mm Bar-%i;15 mm Bar [CH%i+CH%i];10 mm Bar [CH%i+CH%i]",i+1,j+1,2*i,2*i+1,2*j,2*j+1));
			hpl2d[i][j]->Draw();
			myCan->SaveAs("15mmBarTo10mmBar.pdf");
		}
	}
	myCan->SaveAs("15mmBarTo10mmBar.pdf)");
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PlotAllChannels(){
	TH1F *hist[32];
	TH1F *histc[32];
	int bar[32];
	TH1F *m;
	for(int i=0; i<32;i++)
	{
		m = (TH1F*) gROOT->FindObject(Form("htp%i",i));
		if (m) delete m;
	}
	//TPostScript *ps = new TPostScript("AllChannels.ps",112);
	//ps->NewPage();
	for(int i = 0; i < 16; i++)
	{
		c2v->cd(i+1);
		gPad-> SetLogy();
		mppc->Draw(Form("chg[%i]>> htp%i(1000,0,4050)",i,i));
		hist[i] = (TH1F*)gDirectory->Get(Form("htp%i",i));
		histc[i] = (TH1F*)hist[i]->Clone();
		bar[i] = (i+2)/2;
		hist[i]->SetTitle(Form("Single Fibers [CH:%i] on 15 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		hist[i]->SetStats(kFALSE);
		hist[i]->Draw();
	}
	c2v->Update();
	c2v->Print("AllChannels.pdf(");
	for(int i = 16; i < 32; i++)
	{
		c2vt->cd(i-15);
		gPad-> SetLogy();
		mppc->Draw(Form("chg[%i]>> htp%i(1000,0,4050)",i,i));
		hist[i] = (TH1F*)gDirectory->Get(Form("htp%i",i));
		histc[i] = (TH1F*)hist[i]->Clone();
		bar[i] = (i-14)/2;
		hist[i]->SetTitle(Form("Single Fibers [CH%i] on 10 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		hist[i]->SetStats(kFALSE);
		hist[i]->Draw();
	}
	c2vt->Update();
	c2v->Update();
	c2v->Print("AllChannels.pdf)");
	//ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PlotAllFiredChannels()
{
	cout << "Spectrum will be plotted to the channels act in the coincidence!" << endl;
	Mppc *data = filereader.tmCD();
	mppc->Process("FEB_InCoincidenceChannels.C+");
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void SetChannelNumber()
{
	ChannelNumber = fNumberEntry886->GetNumber();
	cout << "Channel number is " << ChannelNumber << endl;
	
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void GetGain()
{
	if(!(fGain->IsOn())) 
	{
		cout << "Gain is not extracted!" << endl;
		ApplyGainFitting = kFALSE;
	}
	if(fGain->IsOn())
	{
		cout << "Gain will be extracted"<< endl;
		ApplyGainFitting = kTRUE;
	}
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Float_t FittingGain2nd(TH1F *h1,TCanvas *myCan)
{
	gStyle->SetOptFit(1);
	Float_t FitGain_gain, FitGain_error;
	TSpectrum *s = new TSpectrum();
	Int_t nfound = s->Search(h1,2,"",0.006); // Find the peaks with 5 sigmas separation from background
	if(nfound>=50) nfound = s->Search(h1,2,"",0.01);
	Int_t npeaks = s->GetNPeaks();                             // Get number of peaks
	Double_t x[100], ex[100];
	Double_t y[100], ey[100];
	Double_t gx[100], gex[100];
	Double_t gy[100], gey[100];
	for (int j=0 ; j<npeaks ; j++)
	{
		x[j] =j+1;
		y[j]=s->GetPositionX()[j];
	}

	//Sort peaks
	Double_t temp;
	int nchanges=0;
	do
	{
		nchanges=0;
		for(int p=0;p<npeaks-1;p++)
		{
			if(y[p]>y[p+1])
			{
				temp=y[p];
				y[p]=y[p+1];
				y[p+1]=temp;
				nchanges++;
			}
		}
	}
	while(nchanges!=0);
	myCan->cd(3);
	int gg = 1;
	for (int g=0 ; g<npeaks ; g++)
	{
		TF1 *gfit = new TF1("gfit","gaus",y[g]-20,y[g]+20);
		h1->Fit(gfit,"QR+");
		if( (abs(y[g] - gfit->GetParameter(1)) ) < 25 )
		{
			gx[gg]  = gg;
			gy[gg]  = gfit->GetParameter(1);
			gey[gg] = gfit->GetParError(1);
			//cout << " --> " << gg << " " << gx[gg] << " " << gy[gg] << " " << gey[gg] << endl;
			gg++;
		}
	}
	myCan->Update();
	TGraphErrors* grpeaks = new TGraphErrors(gg,gx,gy,0,gey);
	for (int p=0; p<gg ; p++)
	{
		grpeaks->SetPoint(p,p,0);
		grpeaks->SetPointError(p,p,0);
	}
	for(int p=1;p<gg;p++)
	{
		grpeaks->SetPoint(p,gx[p],gy[p]);
		grpeaks->SetPointError(p,0,gey[p]);
		//cout << p << " " << gx[p] << " " << gy[p] << " " << gey[p] << endl;
	}
	myCan->cd(4);
	grpeaks->SetTitle("SiPM Gain");
	grpeaks->GetXaxis()->SetTitle("Peak N+i");
	grpeaks->GetYaxis()->SetTitle("ADC Channel");
	grpeaks->SetMarkerColor(4);
	grpeaks->SetMarkerStyle(20);
	grpeaks->SetFillColor(0);
//	TF1 *fit = new TF1("fit","[0] + [1]*x",0.5,3.5);
	TF1 *fit = new TF1("fit","[0] + [1]*x",0.5,4.5);
	fit->SetParName(1,"Gain");
	//fit->SetParName(0, "Pedestal");
	fit->SetParameter(1,80);
	fit->SetParameter(0,200);
	grpeaks->Fit(fit, "QR");
//	gStyle->SetOptFit();
	grpeaks->Draw("AL*");
	myCan->Update();
	Double_t p1 = fit->GetParameter(1);
	Double_t p1error = fit->GetParError(1);
	Double_t p0 = fit->GetParameter(0);
	Double_t p0error = fit->GetParError(0);
	cout<<"Gain "<<p1<<" error "<<p1error<<" ADC Counts/photoelectron"<<endl;
	//cout<<"Pedestal "<<p0<<" error "<<p0error<<endl;
	FitGain_gain=p1;
	if(gg<3) p1error=3.0; //artificially increase error if less than 3 peaks avalable
	FitGain_error=p1error;
	return p1;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Float_t FittingGain(TH1F *h1,TCanvas *myCan, TPostScript *ps)
{
	Float_t FitGain_gain, FitGain_error;
	TSpectrum *s = new TSpectrum();
	Int_t nfound = s->Search(h1,2,"",0.006); // Find the peaks with 5 sigmas separation from background
	if(nfound>=50) nfound = s->Search(h1,2,"",0.01);
	Int_t npeaks = s->GetNPeaks();                             // Get number of peaks
	Double_t x[100], ex[100];
	Double_t y[100], ey[100];
	Double_t gx[100], gex[100];
	Double_t gy[100], gey[100];
	for (int j=0 ; j<npeaks ; j++)
	{
		x[j] =j+1;
		y[j]=s->GetPositionX()[j];
	}

	//Sort peaks
	Double_t temp;
	int nchanges=0;
	do
	{
		nchanges=0;
		for(int p=0;p<npeaks-1;p++)
		{
			if(y[p]>y[p+1])
			{
				temp=y[p];
				y[p]=y[p+1];
				y[p+1]=temp;
				nchanges++;
			}
		}
	}
	while(nchanges!=0);
	int gg = 1;
	for (int g=0 ; g<npeaks ; g++)
	{
		TF1 *gfit = new TF1("gfit","gaus",y[g]-20,y[g]+20);
		h1->Fit(gfit,"QR+");
		if( (abs(y[g] - gfit->GetParameter(1)) ) < 25 )
		{
			gx[gg]  = gg;
			gy[gg]  = gfit->GetParameter(1);
			gey[gg] = gfit->GetParError(1);
			//cout << " --> " << gg << " " << gx[gg] << " " << gy[gg] << " " << gey[gg] << endl;
			gg++;
		}
	}
	myCan->Update();
	ps->NewPage();
	TGraphErrors* grpeaks = new TGraphErrors(gg,gx,gy,0,gey);
	for (int p=0; p<gg ; p++)
	{
		grpeaks->SetPoint(p,p,0);
		grpeaks->SetPointError(p,p,0);
	}
	for(int p=1;p<gg;p++)
	{
		grpeaks->SetPoint(p,gx[p],gy[p]);
		grpeaks->SetPointError(p,0,gey[p]);
		//cout << p << " " << gx[p] << " " << gy[p] << " " << gey[p] << endl;
	}
	grpeaks->SetTitle("SiPM Gain");
	grpeaks->GetXaxis()->SetTitle("Peak N+i");
	grpeaks->GetYaxis()->SetTitle("ADC Channel");
	grpeaks->SetMarkerColor(4);
	grpeaks->SetMarkerStyle(20);
	grpeaks->SetFillColor(0);
//	TF1 *fit = new TF1("fit","[0] + [1]*x",0.5,3.5);
	TF1 *fit = new TF1("fit","[0] + [1]*x",0.5,4.5);
	fit->SetParName(1,"Gain");
	//fit->SetParName(0, "Pedestal");
	fit->SetParameter(1,80);
	fit->SetParameter(0,200);
	grpeaks->Fit(fit, "QR");
//	gStyle->SetOptFit();
	grpeaks->Draw("AL*");
	myCan->Update();
	Double_t p1 = fit->GetParameter(1);
	Double_t p1error = fit->GetParError(1);
	Double_t p0 = fit->GetParameter(0);
	Double_t p0error = fit->GetParError(0);
	cout<<"Gain "<<p1<<" error "<<p1error<<" ADC Counts/photoelectron"<<endl;
	//cout<<"Pedestal "<<p0<<" error "<<p0error<<endl;
	FitGain_gain=p1;
	if(gg<3) p1error=3.0; //artificially increase error if less than 3 peaks avalable
	FitGain_error=p1error;
	return p1;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PlotSingleChannel(){
//	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
//	if (C) delete C;
//	TCanvas *myCan = new TCanvas("Canvas","Quality Check",200,10,700,500);
//	myCan->Divide(4,4);
	if(ApplyGainFitting)
		cout << "WARNING!! Gain will be extracted " << endl;
	TH1F *hist;
	TH1F *histc;
	int bar;
	TH1F *m = (TH1F*) gROOT->FindObject(Form("htp%i",ChannelNumber));
	if (m) delete m;
	c2vs->cd();
	gPad-> SetLogy();
	mppc->Draw(Form("chg[%i]>> htp%i(1000,0,4050)",ChannelNumber,ChannelNumber));
	mppc->Draw(Form("chg[%i]>> htp2%i(1000,250,4050)",ChannelNumber,ChannelNumber));
	hist = (TH1F*)gDirectory->Get(Form("htp%i",ChannelNumber));
	histc = (TH1F*)gDirectory->Get(Form("htp2%i",ChannelNumber));
	bar = (ChannelNumber+2)/2;
	if(ChannelNumber < 16) hist->SetTitle(Form("Single Fibers [CH:%i] on 15 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",ChannelNumber,bar,ChannelNumber));
	if(ChannelNumber > 15) hist->SetTitle(Form("Single Fibers [CH%i] on 10 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",ChannelNumber,bar,ChannelNumber));
	hist->SetStats(kFALSE);
	hist->Draw();
	if(ApplyGainFitting)
	{
		c2vs->Clear();
		c2vs->Modified();
		c2vs->Update();
		c2vs->Divide(2,2);
		c2vs->cd(1);
		hist->Draw();
		c2vs->Update();
		c2vs->cd(2);
		histc->Draw();
		c2vs->Update();
		
		//
		FittingGain2nd(histc,c2vs);
	}
	c2vs->Update();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PlotAllChannels_first(){
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Quality Check",200,10,700,500);
	myCan->Divide(4,4);
	TH1F *hist[32];
	TH1F *histc[32];
	int bar[32];
	TPostScript *ps = new TPostScript("AllChannels.ps",112);
	ps->NewPage();
	for(int i = 0; i < 32; i++)
	{
		
		if(i==16) ps->NewPage();
		if(i<16) myCan->cd(i+1);
		if(i>15) myCan->cd(i-15);
		gPad-> SetLogy();
		mppc->Draw(Form("chg[%i]>> htp%i(1000,0,4050)",i,i));
		hist[i] = (TH1F*)gDirectory->Get(Form("htp%i",i));
		histc[i] = (TH1F*)hist[i]->Clone();
		if(i < 16) bar[i] = (i+2)/2;
		if(i > 15) bar[i] = (i-14)/2;
		cout << "bar number " << bar[i] << endl;
		if(i < 16) hist[i]->SetTitle(Form("Single Fibers [CH:%i] on 15 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		if(i > 15) hist[i]->SetTitle(Form("Single Fibers [CH%i] on 10 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		hist[i]->SetStats(kFALSE);
		hist[i]->Draw();
		if(i==15) myCan->Update();

	}
	myCan->Update();
	myCan->Clear();
	myCan->Divide(1,1);
	for(int i = 0; i < 32; i++)
	{
		myCan->cd(1);
		gPad->SetLogy();
		if(i < 16 ) histc[i]->SetTitle(Form("Single Fibers [CH:%i] on 15 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		if(i > 15 ) histc[i]->SetTitle(Form("Single Fibers [CH%i] on 10 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		histc[i]->SetStats(kFALSE);
		histc[i]->Draw();
		myCan->Update();
		ps->NewPage();
	}
	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PlotGainAllChannels(){
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Quality Check",200,10,700,500);
	TH1F *hist[32];
	Float_t Gain[32];
	int bar[32];
	TPostScript *ps = new TPostScript("AllGain.ps",112);
	ps->NewPage();
	//myCan->Divide(4,4);
	
	for(int i = 0; i < 32; i++)
	{	
		//myCan->cd(i+1);
		gPad-> SetLogy();
		mppc->Draw(Form("chg[%i]>> htpg%i(700,0,4050)",i,i));
		hist[i] = (TH1F*)gDirectory->Get(Form("htpg%i",i));
		if(i < 16) bar[i] = (i+2)/2;
		if(i > 15) bar[i] = (i-14)/2;
		cout << "bar number " << bar[i] << endl;
		if(i < 16) hist[i]->SetTitle(Form("Single Fibers [CH:%i] on 15 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		if(i > 15) hist[i]->SetTitle(Form("Single Fibers [CH%i] on 10 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		//hist[i]->SetStats(kFALSE);
		Gain[i] = FittingGain(hist[i],myCan,ps);
		if(i == 0) myCan->Print("AllGain.pdf(");
		if(i != 0 && i != 31) myCan->Print("AllGain.pdf");
		if(i == 31) myCan->Print("AllGain.pdf)");
	}
	
	ps->Close();
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PlotGainOrChannels(){
	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Quality Check",200,10,700,500);
	TH1F *hist[32];
	Float_t Gain[32];
	int bar[32];
	TPostScript *ps = new TPostScript("AllOrGain.ps",112);
	ps->NewPage();
	//myCan->Divide(4,4);
	for(int i = 0; i < 32; i++)
	{	
		//myCan->cd(i+1);
		gPad-> SetLogy();
		mppc->Draw(Form("chg[%i]>> htpg%i(700,400,4050)",i,i));
		hist[i] = (TH1F*)gDirectory->Get(Form("htpg%i",i));
		if(i < 16) bar[i] = (i+2)/2;
		if(i > 15) bar[i] = (i-14)/2;
		cout << "bar number " << bar[i] << endl;
		if(i < 16) hist[i]->SetTitle(Form("Single Fibers [CH:%i] on 15 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		if(i > 15) hist[i]->SetTitle(Form("Single Fibers [CH%i] on 10 mm Scintillator Bar #%i;Single Fiber on CH%i;Entries",i,bar[i],i));
		//hist[i]->SetStats(kFALSE);
		Gain[i] = FittingGain(hist[i],myCan,ps);
	}
	ps->Close();
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
void TimingDiff(){
	UInt_t time0[500000];
	UInt_t time1[500000];
	Int_t deltaT0, deltaT1;
	//
	Mppc *data = filereader.tmCD();
	Long64_t entries = data->GetInputTree()->GetEntries();
	//
	TH1F *m = (TH1F*) gROOT->FindObject("dT0");
	if (m) delete m;
	TH1F *m2 = (TH1F*) gROOT->FindObject("dT1");
	if (m2) delete m2;

	TH1F *hdT0 = new TH1F("dT0","#DeltaT0",25,-25,25);
	TH1F *hdT1 = new TH1F("dT1","#DeltaT1",25,-25,25);
	hdT0->SetTitle("CAEN FEB (dTS0) ;ns;Entries");
	hdT0->SetLineWidth(2);
	hdT1->SetTitle("CAEN FEB (dTS1);ns;Entries");
	hdT1->SetLineWidth(2);
	TGraph  *gts0;
	TGraph  *gts1;
	TGraph  *gts0ref;
	gts0ref = new TGraph();
	gts1    = new TGraph();
	gts0    = new TGraph();
	//
	for(int i = 0; i < entries; i++)
	{
		data->GetMppc(i);
		time0[i] = data->TS0();
		time1[i] = data->TS1();
		
		//gts0ref->SetPoint(gts0->GetN(),gts0->GetN(),data->TS0Ref()-1e9);
		if(data->TS1()!=0) 
			gts1->SetPoint(gts1->GetN(),gts1->GetN(),data->TS1());
		if(data->TS0()!=0) 
			gts0->SetPoint(gts0->GetN(),gts0->GetN(),data->TS0());
		
		if(i>0)
		{
			deltaT0 = (time0[i] - time0[i-1]);
			hdT0->Fill(deltaT0);
			deltaT1 = (time1[i] - time1[i-1]);
			hdT1->Fill(deltaT1);
		}
	}
	myCant->cd(1);
	hdT0->Draw();
	myCant->Update();
	myCant->cd(2);
	hdT1->Draw();
	myCant->Update();
	myCant->Modified();
	myCant->cd(3);
	gts0ref->Draw("AL");
	gts0ref->GetHistogram()->GetXaxis()->SetTitle("Event number");
	gts0ref->GetHistogram()->GetYaxis()->SetTitle("TS0 period deviation from 1s, ns");
	myCant->Update();
	myCant->Modified();
	myCant->cd(4);
	gts1->SetLineColor(kRed);
	gts1->Draw("AL");
	myCant->Update();
	myCant->Modified();
	gts1->GetHistogram()->GetXaxis()->SetTitle("Event number");
	gts1->GetHistogram()->GetYaxis()->SetTitle("TS0, TS1 (red), ns");
	gts0->Draw("sameL"); 
	myCant->Update();
	myCant->Modified();
}
*/
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void CheckCoincidenceSingleNew(int ientry){
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLabelOffset(0.007,"X");
	gStyle->SetLabelOffset(0.007,"Y");
	TPostScript *ps = new TPostScript(Form("Coincidence_entry%i.ps",ientry),112);
	TH1F *histos = new TH1F("histos",Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);;
	TH1F *h_CAENCoin = new TH1F("h_CAENCoin",Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);;
	TH1F *h_PosCoinBottom = new TH1F("h_PosCoinBottom",Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);;
	TH1F *h_PosCoinTop = new TH1F("h_PosCoinTop",Form("Cosmic Muon Event %i;Channel Number;ADC",ientry),32,-0.5,31.5);;
	Int_t ScintillatorX;
	Int_t ScintillatorY;
	
	TH2F *hCoinPlane = new TH2F("hCoinPlane","Cosmic Muon Event in Coincidence;15 mm Scintillator Bar;10 mm Scintillator Bar",20,0,8,20,0,8);

	TCanvas *C = (TCanvas*) gROOT->FindObject("Canvas");
	if (C) delete C;
	TCanvas *myCan = new TCanvas("Canvas","Quality Check",200,10,700,500);

	Mppc *data = filereader.tmCD();
	UInt_t BarX, BarY;
	Int_t num;
	Int_t InCoincidence;
	std::vector<int> v_InCoincidence;
	UInt_t PossibleCoincidence;
	UInt_t PossibleCoincidence2;
	UInt_t MaxAdcInPosCoincidence;
	UInt_t MinAdcInPosCoincidence;
	UInt_t MaxAdcInPosCoincidence2;
	UInt_t MinAdcInPosCoincidence2;

	Double_t ratioBottom[8];
	Double_t weightedBottom[8];
	Double_t ratioTop[8];
	Double_t weightedTop[8];
	TBox *LayerX, *LayerY;
	//
	data->GetMppc(ientry);
	for (int i=0;i<8;i++)
	{
		if (data->ADC(2*i) == 0 || data->ADC(2*i+1) ==0)
		{
			ratioBottom[i]    = 0;
			weightedBottom[i] = 0;
		}else{
			ratioBottom[i] = (Double_t)data->ADC(2*i)/data->ADC(2*i+1);
			weightedBottom[i] = (ratioBottom[i]+1)*(data->ADC(2*i+1)+data->ADC(2*i));
			if(ratioBottom[i]>1.) {
				ratioBottom[i] = (Double_t)data->ADC(2*i+1)/data->ADC(2*i);
				weightedBottom[i] = (ratioBottom[i]+1)*(data->ADC(2*i)+data->ADC(2*i+1));
			}
		}
		if(DEBUG1) 
			cout << "Bottom ----> ratio for CH#" << 2*i << "-CH#" << 2*i+1 << " " << ratioBottom[i] << " " << weightedBottom[i] << " " << endl;
		//
		if (data->ADC(2*i+16) == 0 || data->ADC(2*i+1+16) ==0)
		{
			ratioTop[i]    = 0;
			weightedTop[i] = 0;
		}else{
			ratioTop[i] = (Double_t)data->ADC(2*i+16)/data->ADC(2*i+1+16);
			weightedTop[i] = (ratioTop[i]+1)*(data->ADC(2*i+1+16)+data->ADC(2*i+16));
		if(ratioTop[i]>1.) {
				ratioTop[i] = (Double_t)data->ADC(2*i+1+16)/data->ADC(2*i+16);
				weightedTop[i] = (ratioTop[i]+1)*(data->ADC(2*i+16)+data->ADC(2*i+1+16));
			}
		}
		if(DEBUG1) 
			cout << "Top ----> ratio for CH#" << 2*i+16 << "-CH#" << 2*i+1+16 << " " << ratioTop[i] << " " << weightedTop[i] << " " << endl;
	}
	PossibleCoincidence = TMath::LocMax(8,weightedBottom);
	MaxAdcInPosCoincidence = TMath::Max(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
	MinAdcInPosCoincidence = TMath::Min(data->ADC(2*PossibleCoincidence),data->ADC(2*PossibleCoincidence+1));
	if(DEBUG1) 
	cout << ientry 
		 << " locate max. " 
		 << PossibleCoincidence 
		 << " " << weightedBottom[PossibleCoincidence] 
		 << " CH#" << 2*PossibleCoincidence 
		 << "-CH#" << 2*PossibleCoincidence+1 
		 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence)
		 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence+1)
		 << " maxADC " << MaxAdcInPosCoincidence
		 << " minADC " << MinAdcInPosCoincidence
		 << endl;
	
	ScintillatorY = PossibleCoincidence;
	
	Long64_t TwoChannelsInPosCoincidence_Plane1 = (2*PossibleCoincidence+2*PossibleCoincidence+1)/2;
	if(DEBUG1) cout << "TwoChannelsInPosCoincidence_Plane1 " << TwoChannelsInPosCoincidence_Plane1 << endl;
	//
	PossibleCoincidence2 = TMath::LocMax(8,weightedTop);
	MaxAdcInPosCoincidence2 = TMath::Max(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
	MinAdcInPosCoincidence2 = TMath::Min(data->ADC(2*PossibleCoincidence2+16),data->ADC(2*PossibleCoincidence2+1+16));
	if(DEBUG1) 
	cout << ientry 
		 << " locate max2. " 
		 << PossibleCoincidence2 
		 << " " << weightedTop[PossibleCoincidence2] 
		 << " CH#" << 2*PossibleCoincidence2+16
		 << "-CH#" << 2*PossibleCoincidence2+1+16
		 << " ADC[CH1st] " << data->ADC(2*PossibleCoincidence2+16)
		 << " ADC[CH2nd] " << data->ADC(2*PossibleCoincidence2+1+16)
		 << " maxADC " << MaxAdcInPosCoincidence2
		 << " minADC " << MinAdcInPosCoincidence2
		 << endl;

	hCoinPlane->Fill(PossibleCoincidence2+1,PossibleCoincidence+1);
		
	ScintillatorX = PossibleCoincidence2;
	Long64_t TwoChannelsInPosCoincidence_Plane2 = (2*PossibleCoincidence2+16+2*PossibleCoincidence2+1+16)/2;
	if(DEBUG1) cout << "TwoChannelsInPosCoincidence_Plane2 " << TwoChannelsInPosCoincidence_Plane2 << endl;
	for (int i = 0; i < 32; i++ )
	{
		cout << i << " " << data->ADC(i) << " " << data->Coincidence() << " " << endl;
		//if(i==31) InCoincidence = ConvertToBinary(ientry,data->Coincidence(),16);
		if(i==31) v_InCoincidence = ConvertToBinary2nd(ientry,data->Coincidence(),16);
		cout << endl;
		histos->SetBinContent(i+1,data->ADC(i));
	}
	cout << "UMUTTTTT --> " << v_InCoincidence.size() << endl;
	
	//if(DEBUG) cout << "Coincidence " << InCoincidence << " " << PossibleCoincidence << " " << PossibleCoincidence2 << endl;
	//if(DEBUG) cout << "CAEN coincidence pl " << InCoincidence << " " << 2*InCoincidence << " " << 2*InCoincidence+1 << endl;
	//h_CAENCoin->SetBinContent(2*InCoincidence-1,data->ADC(2*InCoincidence-2));
	//h_CAENCoin->SetBinContent(2*InCoincidence,data->ADC(2*InCoincidence-1));

	for(int i=0; i<v_InCoincidence.size();i++)
	{
		if(DEBUG1) cout << "Coincidence " << v_InCoincidence[i] << " " << PossibleCoincidence << " " << PossibleCoincidence2 << endl;
		if(DEBUG1) cout << "CAEN coincidence pl " << v_InCoincidence[i] << " " << 2*v_InCoincidence[i] << " " << 2*v_InCoincidence[i]+1 << endl;
		h_CAENCoin->SetBinContent(2*v_InCoincidence[i]-1,data->ADC(2*v_InCoincidence[i]-2));
		h_CAENCoin->SetBinContent(2*v_InCoincidence[i],data->ADC(2*v_InCoincidence[i]-1));
	}
	
	
	h_PosCoinBottom->SetBinContent(TwoChannelsInPosCoincidence_Plane1+1,data->ADC(TwoChannelsInPosCoincidence_Plane1));
	h_PosCoinBottom->SetBinContent(TwoChannelsInPosCoincidence_Plane1+2,data->ADC(TwoChannelsInPosCoincidence_Plane1+1));
	h_PosCoinTop->SetBinContent(TwoChannelsInPosCoincidence_Plane2+1,data->ADC(TwoChannelsInPosCoincidence_Plane2));
	h_PosCoinTop->SetBinContent(TwoChannelsInPosCoincidence_Plane2+2,data->ADC(TwoChannelsInPosCoincidence_Plane2+1));

	histos->SetMaximum(5000);
	histos->Draw();
	histos->GetYaxis()->SetTitleOffset(1.3);
	h_CAENCoin->SetFillColor(kRed);
	h_CAENCoin->SetFillStyle(1001);
	h_CAENCoin->Draw("same");
	h_PosCoinBottom->SetFillColor(kGreen);
	h_PosCoinBottom->SetFillStyle(3012);
	h_PosCoinBottom->Draw("same");
	h_PosCoinTop->SetFillColor(kGreen);
	h_PosCoinTop->SetFillStyle(3022);
	h_PosCoinTop->Draw("same");
	TLine *line = new TLine(15.5,0,15.5,5000);
	line->SetLineColor(kBlue);
	line->SetLineWidth(2);
	line->Draw("same");
	TText *t = new TText(7.5,4500,"Bottom Layer (15mm Bars)");
	t->SetTextAlign(22);
	t->SetTextColor(kBlue);
	t->SetTextFont(43);
	t->SetTextSize(14);
	t->Draw("same");
	TText *t2 = new TText(23,4500,"Top Layer (10mm Bars)");
	t2->SetTextAlign(22);
	t2->SetTextColor(kBlue);
	t2->SetTextFont(43);
	t2->SetTextSize(14);
	t2->Draw("same");
	myCan->Update();
	//ps->NewPage();
/*
	TH1 *frame = myCan->DrawFrame(0,0,1850,1850);
	frame->GetXaxis()->SetBinLabel(80,"Scin8");
	frame->GetXaxis()->SetBinLabel(200,"Scin7");
	frame->GetXaxis()->SetBinLabel(320,"Scin6");
	frame->GetXaxis()->SetBinLabel(440,"Scin5");
	frame->GetXaxis()->SetBinLabel(560,"Scin4");
	frame->GetXaxis()->SetBinLabel(680,"Scin3");
	frame->GetXaxis()->SetBinLabel(800,"Scin2");
	frame->GetXaxis()->SetBinLabel(920,"Scin1");
	frame->GetXaxis()->SetLabelSize(0.04);

	frame->GetYaxis()->Set(8,frame->GetMinimum(),frame->GetMaximum());
	frame->GetYaxis()->SetBinLabel(8,"Scin16");
	frame->GetYaxis()->SetBinLabel(7,"Scin15");
	frame->GetYaxis()->SetBinLabel(6,"Scin14");
	frame->GetYaxis()->SetBinLabel(5,"Scin13");
	frame->GetYaxis()->SetBinLabel(4,"Scin12");
	frame->GetYaxis()->SetBinLabel(3,"Scin11");
	frame->GetYaxis()->SetBinLabel(2,"Scin10");
	frame->GetYaxis()->SetBinLabel(1,"Scin9");
	frame->GetYaxis()->SetLabelSize(0.04);

	frame->GetXaxis()->SetTitleOffset(1.6);
	frame->GetXaxis()->SetTitle("Bottom Layer (15 mm)");
	frame->GetYaxis()->SetTitleOffset(1.6);
	frame->GetYaxis()->SetTitle("Top Layer (10 mm)");
	frame->SetTitle("CRT Module (Scintillators in Coincidence)");
	//
	LayerY = new TBox(1850-230*ScintillatorY,0,1850-230*(ScintillatorY+1),1850);
	LayerY->SetFillColor(kGray+1);
	LayerY->SetFillStyle(3006);
	LayerY->SetLineColor(kBlue);
	LayerY->SetLineWidth(2);
	LayerY->Draw();
	LayerX = new TBox(0,230*ScintillatorX,1850,230*ScintillatorX+230);
	LayerX->SetFillColor(kGray+1);
	LayerX->SetFillStyle(3007);
	LayerX->SetLineColor(kRed);
	LayerX->SetLineWidth(2);
	LayerX->Draw("same");
	myCan->Update();
*/
	ps->Close();
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Float_t FittingPedestal(TH1F *hp)
{
	Float_t mean, sigma;
	Float_t maxadc=hp->GetXaxis()->GetBinCenter(hp->GetMaximumBin());
	TF1 *gfit = new TF1("gfit","gaus",maxadc-10,maxadc+10);
	hp->Fit(gfit,"R+");
	mean=gfit->GetParameter(1);
	sigma=gfit->GetParameter(2);
	hp->Fit(gfit,"","",mean-2*sigma,mean+2*sigma);
	cout<<"Pedestal "<< mean <<" error "<< sigma <<endl;
	return mean;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void SetEventNumber()
{
	EventNumber = fNumberEntry75->GetNumber();
	cout << "Event Number is " << EventNumber << endl;
	
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void CurrentCRTEvent()
{
	CurrentEventNumber = EventNumber;
	ShowEventOnCRT_CAEN(EventNumber);
	CheckCoincidenceSingle(EventNumber);
	cout<< "Current Event Number is " << EventNumber << endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void NextCRTEvent()
{
	CurrentEventNumber = EventNumber;
	EventNumber = CurrentEventNumber + 1;
	ShowEventOnCRT_CAEN(EventNumber);
	CheckCoincidenceSingle(EventNumber);
	cout<< "Showing Event Number is " << EventNumber << endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PreviousCRTEvent()
{
	CurrentEventNumber = EventNumber;
	EventNumber = CurrentEventNumber - 1;
	ShowEventOnCRT_CAEN(EventNumber);
	CheckCoincidenceSingle(EventNumber);
	cout<< "Showing Event Number is " << EventNumber << endl;
	
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void NextScinEvent()
{
	CurrentEventNumber = EventNumber;
	EventNumber = CurrentEventNumber + 1;
	CheckCoincidenceSingle(EventNumber);
	cout<< "Showing Event Number is " << EventNumber << endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void PreviousScinEvent()
{
	CurrentEventNumber = EventNumber;
	EventNumber = CurrentEventNumber - 1;
	CheckCoincidenceSingle(EventNumber);
	cout<< "Showing Event Number is " << EventNumber << endl;
	
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void Exit()
{
	cout << "Exit application..." << endl;
	gROOT->Reset();
	gApplication->Terminate(0);
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vector<int> noise_v;

void SetNoisyEvent()
{
	Int_t ev;
	FILE *fp;
	char buf[512];
	if ( (fp = fopen( "NoiseList.txt", "rt" )) == NULL) {
	        cout << endl << "File does not exist!!!" << endl;
	        exit(0);
	}	
	while (fgets( buf, sizeof(buf), fp ) != 0) {
		sscanf(buf,"%d",&ev);
		noise_v.push_back(ev);
	}
	fclose(fp);
	cout << "Number of noisy event located in the file:  " << noise_v.size() << endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Int_t NoisyEventNumber;
Int_t CurrentNoisyEventNumber;
void ShowNoisyEvent()
{
	CurrentNoisyEventNumber = 0;
	CheckCoincidenceSingle(noise_v[CurrentNoisyEventNumber]);
	cout<< "Showing Event Number is " << CurrentNoisyEventNumber << endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ShowNextNoisyEvent()
{
	NoisyEventNumber = CurrentNoisyEventNumber;
	CurrentNoisyEventNumber = NoisyEventNumber + 1;
	CheckCoincidenceSingle(noise_v[CurrentNoisyEventNumber]);
	cout<< "Showing Event Number is " << noise_v[CurrentNoisyEventNumber] << endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void ShowPreviousNoisyEvent()
{
	NoisyEventNumber = CurrentNoisyEventNumber;
	CurrentNoisyEventNumber = NoisyEventNumber - 1;
	CheckCoincidenceSingle(noise_v[CurrentNoisyEventNumber]);
	cout<< "Showing Event Number is " << noise_v[CurrentNoisyEventNumber] << endl;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void FEBSetting()
{
	TVectorD *tB  = (TVectorD*) _file0->Get("infoBias");
	TVectorD *tP  = (TVectorD*) _file0->Get("infoPreAmp");
	TNamed *tN = (TNamed*) _file0->Get("BarcodeInfo");
	TNamed *tN2 = (TNamed*) _file0->Get("StartT");
	TNamed *tN3 = (TNamed*) _file0->Get("StopT");
	TNamed *tN4 = (TNamed*) _file0->Get("DACInfo");
	TNamed *tN5 = (TNamed*) _file0->Get("TRateInfo");
	if(tB!=NULL||tP!=NULL||tN!=NULL)
	{
		myCants->cd();
		myCants->Clear();
		TPaveText *pV = new TPaveText(0.01,0.01,0.955,0.955,"NDC");
		pV->SetTextFont(32);
		pV->SetTextColor(4);
		pV->SetTextAlign(12);
		pV->AddText("");
		pV->AddText(Form("File Name %s",_file0->GetName()));
		
		pV->AddText(Form("Barcode number %s",tN->GetName()));
		
		pV->AddText(Form("DAQ starting time %s",tN2->GetName()));
		pV->AddText(Form("DAQ stopping time %s",tN3->GetName()));
		pV->AddText(Form("Threshould %s",tN4->GetName()));
		pV->AddText(Form("Trigger Rate %s",tN5->GetName()));
		
		pV->AddText("The measurement has been performed with following configurations:");
		for(Int_t i=0; i<tB->GetNrows(); i++)
		{
			pV->AddText(Form(" CH#%i          %.0f         %.0f",i, (*tP)[i],(*tB)[i]));
		}
		pV->SetLabel("CRT Module");
		pV->Draw();
		myCants->Modified();
		myCants->Update();
	}
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void DisplayCRT()
{
	// main frame
	TGMainFrame *fMainFrame1314 = new TGMainFrame(gClient->GetRoot(),10,10,kMainFrame | kVerticalFrame);
	fMainFrame1314->SetName("fMainFrame1314");
	fMainFrame1314->SetLayoutBroken(kTRUE);
	fMainFrame1314->SetCleanup(kDeepCleanup);
	// composite frame
	TGCompositeFrame *fMainFrame1560 = new TGCompositeFrame(fMainFrame1314,1329,789,kVerticalFrame);
	fMainFrame1560->SetName("fMainFrame1560");
	fMainFrame1560->SetLayoutBroken(kTRUE);

	// composite frame
	TGCompositeFrame *fMainFrame1241 = new TGCompositeFrame(fMainFrame1560,1329,789,kVerticalFrame);
	fMainFrame1241->SetName("fMainFrame1241");

	// vertical frame
	TGVerticalFrame *fVerticalFrame734 = new TGVerticalFrame(fMainFrame1241,1327,787,kVerticalFrame);
	fVerticalFrame734->SetName("fVerticalFrame734");

	// status bar
	Int_t parts[] = {15, 15, 15, 15, 15, 25};
	fStatusBar739 = new TGStatusBar(fVerticalFrame734,1327,18);
	fStatusBar739->SetName("fStatusBar739");
	fStatusBar739->SetParts(parts, 6);
	fVerticalFrame734->AddFrame(fStatusBar739, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));

	// horizontal frame
	TGHorizontalFrame *fHorizontalFrame768 = new TGHorizontalFrame(fVerticalFrame734,1350,765,kHorizontalFrame);
	fHorizontalFrame768->SetName("fHorizontalFrame768");

	// "DAQ FEB controls" group frame
	TGGroupFrame *fGroupFrame679 = new TGGroupFrame(fHorizontalFrame768,"CRT Module Monitor");

	fGroupFrame679->SetLayoutManager(new TGVerticalLayout(fGroupFrame679));
	fGroupFrame679->Resize(155,761);
	fHorizontalFrame768->AddFrame(fGroupFrame679, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY,2,2,2,2));

	///////////////////////////////////////////////////////////////////////////////////

	fNumberEntry75 = new TGNumberEntry(fGroupFrame679, (Double_t) 0,6,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 1,(TGNumberFormat::ELimit) 2,0,10000);
	fGroupFrame679->AddFrame(fNumberEntry75, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,32,2));
	TGTextButton *fTextButton93 = new TGTextButton(fGroupFrame679,"Event Number");
	fTextButton93->SetTextJustify(36);
	fTextButton93->SetMargins(0,0,0,0);
	fTextButton93->SetWrapLength(-1);
	fTextButton93->Resize(123,22);
	fTextButton93->SetCommand("SetEventNumber()");
	fGroupFrame679->AddFrame(fTextButton93, new TGLayoutHints(kLHintsLeft | kLHintsTop,0,0,2,2));
/*
	TGTextButton *fTextButton788 = new TGTextButton(fGroupFrame679,"Plot Waveforms");
	fTextButton788->SetTextJustify(36);
	fTextButton788->SetMargins(0,0,0,0);
	fTextButton788->SetWrapLength(-1);
	fTextButton788->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton788, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton788->SetCommand("Waveforms()");
*/
	fLabel79 = new TGLabel(fGroupFrame679,"::::::::::: Show on CRT Module :::::::::");
	fLabel79->SetTextJustify(36);
	fLabel79->SetMargins(0,0,0,0);
	fLabel79->SetWrapLength(-1);
	fGroupFrame679->AddFrame(fLabel79, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,10,2));

	TGTextButton *fTextButton787s = new TGTextButton(fGroupFrame679,"Show Event");
	fTextButton787s->SetTextJustify(36);
	fTextButton787s->SetMargins(0,0,0,0);
	fTextButton787s->SetWrapLength(-1);
	fTextButton787s->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton787s, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton787s->SetCommand("CurrentCRTEvent()");

	TGTextButton *fTextButton787 = new TGTextButton(fGroupFrame679,"Show Next Event");
	fTextButton787->SetTextJustify(36);
	fTextButton787->SetMargins(0,0,0,0);
	fTextButton787->SetWrapLength(-1);
	fTextButton787->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton787, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton787->SetCommand("NextCRTEvent()");


	TGTextButton *fTextButton787x = new TGTextButton(fGroupFrame679,"Show Previous Event");
	fTextButton787x->SetTextJustify(36);
	fTextButton787x->SetMargins(0,0,0,0);
	fTextButton787x->SetWrapLength(-1);
	fTextButton787x->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton787x, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton787x->SetCommand("PreviousCRTEvent()");


	TGTextButton *fTextButton787xx = new TGTextButton(fGroupFrame679,"Muon Hits on CRT");
	fTextButton787xx->SetTextJustify(36);
	fTextButton787xx->SetMargins(0,0,0,0);
	fTextButton787xx->SetWrapLength(-1);
	fTextButton787xx->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton787xx, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,10,2));
	fTextButton787xx->SetCommand("CRTMuonHits()");



/*
	fLabel79s = new TGLabel(fGroupFrame679,"::::::::::: Show FEB Channels :::::::::");
	fLabel79s->SetTextJustify(36);
	fLabel79s->SetMargins(0,0,0,0);
	fLabel79s->SetWrapLength(-1);
	fGroupFrame679->AddFrame(fLabel79s, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,20,2));

	TGTextButton *fTextButton787s = new TGTextButton(fGroupFrame679,"Show Next Event");
	fTextButton787s->SetTextJustify(36);
	fTextButton787s->SetMargins(0,0,0,0);
	fTextButton787s->SetWrapLength(-1);
	fTextButton787s->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton787s, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton787s->SetCommand("NextScinEvent()");


	TGTextButton *fTextButton787s = new TGTextButton(fGroupFrame679,"Show Previous Event");
	fTextButton787s->SetTextJustify(36);
	fTextButton787s->SetMargins(0,0,0,0);
	fTextButton787s->SetWrapLength(-1);
	fTextButton787s->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton787s, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton787s->SetCommand("PreviousScinEvent()");
*/
	fLabel79s = new TGLabel(fGroupFrame679,"::::::::::: FEB Behaviour :::::::::");
	fLabel79s->SetTextJustify(36);
	fLabel79s->SetMargins(0,0,0,0);
	fLabel79s->SetWrapLength(-1);
	fGroupFrame679->AddFrame(fLabel79s, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,10,2));

	TGTextButton *fTextButton788r = new TGTextButton(fGroupFrame679,"Performance");
	fTextButton788r->SetTextJustify(36);
	fTextButton788r->SetMargins(0,0,0,0);
	fTextButton788r->SetWrapLength(-1);
	fTextButton788r->Resize(123,22);
	fTextButton788r->SetCommand("FEB_Behaviour_WithPlots()");
	fGroupFrame679->AddFrame(fTextButton788r, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fLabel79s = new TGLabel(fGroupFrame679,"...................................");
	fLabel79s->SetTextJustify(36);
	fLabel79s->SetMargins(0,0,0,0);
	fLabel79s->SetWrapLength(-1);
	fGroupFrame679->AddFrame(fLabel79s, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,2,2));

	TGTextButton *fTextButton787n = new TGTextButton(fGroupFrame679,"Set Noisy Event");
	fTextButton787n->SetTextJustify(36);
	fTextButton787n->SetMargins(0,0,0,0);
	fTextButton787n->SetWrapLength(-1);
	fTextButton787n->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton787n, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton787n->SetCommand("SetNoisyEvent()");


	TGTextButton *fTextButton788n = new TGTextButton(fGroupFrame679,"Show Noisy Event");
	fTextButton788n->SetTextJustify(36);
	fTextButton788n->SetMargins(0,0,0,0);
	fTextButton788n->SetWrapLength(-1);
	fTextButton788n->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton788n, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton788n->SetCommand("ShowNoisyEvent()");

	TGTextButton *fTextButton789n = new TGTextButton(fGroupFrame679,"Show Next Noisy Event");
	fTextButton789n->SetTextJustify(36);
	fTextButton789n->SetMargins(0,0,0,0);
	fTextButton789n->SetWrapLength(-1);
	fTextButton789n->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton789n, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton789n->SetCommand("ShowNextNoisyEvent()");

	TGTextButton *fTextButton790n = new TGTextButton(fGroupFrame679,"Show Previous Noisy Event");
	fTextButton790n->SetTextJustify(36);
	fTextButton790n->SetMargins(0,0,0,0);
	fTextButton790n->SetWrapLength(-1);
	fTextButton790n->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton790n, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton790n->SetCommand("ShowPreviousNoisyEvent()");

	fLabel79s = new TGLabel(fGroupFrame679,"::::::::::: Light Yield Check :::::::::");
	fLabel79s->SetTextJustify(36);
	fLabel79s->SetMargins(0,0,0,0);
	fLabel79s->SetWrapLength(-1);
	fGroupFrame679->AddFrame(fLabel79s, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,10,2));

	TGTextButton *fTextButton788v = new TGTextButton(fGroupFrame679,"All Channels");
	fTextButton788v->SetTextJustify(36);
	fTextButton788v->SetMargins(0,0,0,0);
	fTextButton788v->SetWrapLength(-1);
	fTextButton788v->Resize(123,22);
	fTextButton788v->SetCommand("PlotAllChannels()");
	fGroupFrame679->AddFrame(fTextButton788v, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));

	TGTextButton *fTextButton788vc = new TGTextButton(fGroupFrame679,"All Channels in Coincidence");
	fTextButton788vc->SetTextJustify(36);
	fTextButton788vc->SetMargins(0,0,0,0);
	fTextButton788vc->SetWrapLength(-1);
	fTextButton788vc->Resize(123,22);
	fTextButton788vc->SetCommand("PlotAllFiredChannels()");
	fGroupFrame679->AddFrame(fTextButton788vc, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));


	TGLabel *fLabel764 = new TGLabel(fGroupFrame679,"Set Channel #");
	fLabel764->SetTextJustify(36);
	fLabel764->SetMargins(0,0,0,0);
	fLabel764->SetWrapLength(-1);
	fGroupFrame679->AddFrame(fLabel764, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX,2,2,10,2));
	fNumberEntry886 = new TGNumberEntry(fGroupFrame679, (Double_t) 4,3,-1,(TGNumberFormat::EStyle) 5,(TGNumberFormat::EAttribute) 1,(TGNumberFormat::ELimit) 2,0,31);
	fGroupFrame679->AddFrame(fNumberEntry886, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX ));
	fNumberEntry886->Connect("ValueSet(Long_t)", 0, 0,  "SetChannelNumber()");


	fGain = new TGCheckButton(fGroupFrame679,"Extract Gain");
	fGain->SetTextJustify(36);
	fGain->SetMargins(0,0,0,0);
	fGain->SetWrapLength(-1);
	fGain->SetCommand("GetGain()");
	fGroupFrame679->AddFrame(fGain, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,2));

	TGTextButton *fTextButton78911 = new TGTextButton(fGroupFrame679,"Single Channel");
	fTextButton78911->SetTextJustify(36);
	fTextButton78911->SetMargins(0,0,0,0);
	fTextButton78911->SetWrapLength(-1);
	fTextButton78911->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton78911, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton78911->SetCommand("PlotSingleChannel()");



	TGTextButton *fTextButton788f = new TGTextButton(fGroupFrame679,"Two fibes in a ScinBar");
	fTextButton788f->SetTextJustify(36);
	fTextButton788f->SetMargins(0,0,0,0);
	fTextButton788f->SetWrapLength(-1);
	fTextButton788f->Resize(123,22);
	fTextButton788f->SetCommand("TwoFibersInBar()");
	fGroupFrame679->AddFrame(fTextButton788f, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));

	TGTextButton *fTextButton788fp = new TGTextButton(fGroupFrame679,"ScinBars in Two Layers");
	fTextButton788fp->SetTextJustify(36);
	fTextButton788fp->SetMargins(0,0,0,0);
	fTextButton788fp->SetWrapLength(-1);
	fTextButton788fp->Resize(123,22);
	fTextButton788fp->SetCommand("BarsInTwoPlanes()");
	fGroupFrame679->AddFrame(fTextButton788fp, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));

	/*	TGTextButton *fTextButton788tt = new TGTextButton(fGroupFrame679,"Timing");
	fTextButton788tt->SetTextJustify(36);
	fTextButton788tt->SetMargins(0,0,0,0);
	fTextButton788tt->SetWrapLength(-1);
	fTextButton788tt->Resize(123,22);
	fTextButton788tt->SetCommand("TimingDiff()");
	fGroupFrame679->AddFrame(fTextButton788tt, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	*/
	TGTextButton *fTextButton787ss = new TGTextButton(fGroupFrame679,"FEB Configurations");
	fTextButton787ss->SetTextJustify(36);
	fTextButton787ss->SetMargins(0,0,0,0);
	fTextButton787ss->SetWrapLength(-1);
	fTextButton787ss->Resize(123,22);
	fGroupFrame679->AddFrame(fTextButton787ss, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,2,2));
	fTextButton787ss->SetCommand("FEBSetting()");

	TGTextButton *fTextButton788 = new TGTextButton(fGroupFrame679,"&Exit");
	fTextButton788->SetTextJustify(36);
	fTextButton788->SetMargins(0,0,0,0);
	fTextButton788->SetWrapLength(-1);
	fTextButton788->Resize(123,22);
	fTextButton788->SetCommand("Exit()");
	fGroupFrame679->AddFrame(fTextButton788, new TGLayoutHints(kLHintsLeft| kLHintsCenterX  | kLHintsTop | kLHintsExpandX,0,0,10,2));


	//////////////////////////////////////////////////////////////////////
	// tab widget
	TGTab *fTab683 = new TGTab(fHorizontalFrame768,1187,761);
	// container of "All histos"
	TGCompositeFrame *fCompositeFrame720;
	fCompositeFrame720 = fTab683->AddTab("Event on CRT");
	fCompositeFrame720->SetLayoutManager(new TGVerticalLayout(fCompositeFrame720));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas721 = new TRootEmbeddedCanvas(0,fCompositeFrame720,1179,732);
	Int_t wfRootEmbeddedCanvas721 = fRootEmbeddedCanvas721->GetCanvasWindowId();
	c = new TCanvas("c", 10, 10, wfRootEmbeddedCanvas721); 	
	fRootEmbeddedCanvas721->AdoptCanvas(c);
	fCompositeFrame720->AddFrame(fRootEmbeddedCanvas721, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));
	//
	// container of "One channel"
	TGCompositeFrame *fCompositeFrame735;
	fCompositeFrame735 = fTab683->AddTab("Coincidence");
	fCompositeFrame735->SetLayoutManager(new TGVerticalLayout(fCompositeFrame735));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas736 = new TRootEmbeddedCanvas(0,fCompositeFrame735,1179,732);
	Int_t wfRootEmbeddedCanvas736 = fRootEmbeddedCanvas736->GetCanvasWindowId();
	c1 = new TCanvas("c1", 10, 10, wfRootEmbeddedCanvas736);
	fRootEmbeddedCanvas736->AdoptCanvas(c1);
	fCompositeFrame735->AddFrame(fRootEmbeddedCanvas736, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));


	// container of "TChannel profile"
	TGCompositeFrame *fCompositeFrame7352;
	fCompositeFrame7352 = fTab683->AddTab("Performance");
	fCompositeFrame7352->SetLayoutManager(new TGVerticalLayout(fCompositeFrame7352));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas7362 = new TRootEmbeddedCanvas(0,fCompositeFrame7352,1179,732);
	Int_t wfRootEmbeddedCanvas7362 = fRootEmbeddedCanvas7362->GetCanvasWindowId();
	myCan = new TCanvas("myCan", 10, 10, wfRootEmbeddedCanvas7362);
	myCan->Divide(2,2);
	fRootEmbeddedCanvas7362->AdoptCanvas(myCan);
	fCompositeFrame7352->AddFrame(fRootEmbeddedCanvas7362, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));
//	c2->SetGridx(1);
//	c2->SetGridy(1);

	// container of "TChannel profile"
	TGCompositeFrame *fCompositeFrame7352x;
	fCompositeFrame7352x = fTab683->AddTab("Flagged PL");
	fCompositeFrame7352x->SetLayoutManager(new TGVerticalLayout(fCompositeFrame7352x));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas7362x = new TRootEmbeddedCanvas(0,fCompositeFrame7352x,1179,732);
	Int_t wfRootEmbeddedCanvas7362x = fRootEmbeddedCanvas7362x->GetCanvasWindowId();
	myCanx = new TCanvas("myCanx", 10, 10, wfRootEmbeddedCanvas7362x);
	myCanx->Divide(1,1);
	fRootEmbeddedCanvas7362x->AdoptCanvas(myCanx);
	fCompositeFrame7352x->AddFrame(fRootEmbeddedCanvas7362x, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));
//	c2->SetGridx(1);
//	c2->SetGridy(1);

	// container of "One channel"
	TGCompositeFrame *fCompositeFrame7353v;
	fCompositeFrame7353v = fTab683->AddTab("All CH [15mm]");
	fCompositeFrame7353v->SetLayoutManager(new TGVerticalLayout(fCompositeFrame7353v));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas7363v = new TRootEmbeddedCanvas(0,fCompositeFrame7353v,1179,732);
	Int_t wfRootEmbeddedCanvas7363v = fRootEmbeddedCanvas7363v->GetCanvasWindowId();
	c2v = new TCanvas("c2v", 10, 10, wfRootEmbeddedCanvas7363v);
	c2v->Divide(4,4);
	fRootEmbeddedCanvas7363v->AdoptCanvas(c2v);
	fCompositeFrame7353v->AddFrame(fRootEmbeddedCanvas7363v, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

	// container of "One channel"
	TGCompositeFrame *fCompositeFrame7353vt;
	fCompositeFrame7353vt = fTab683->AddTab("All CH [10mm]");
	fCompositeFrame7353vt->SetLayoutManager(new TGVerticalLayout(fCompositeFrame7353vt));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas7363vt = new TRootEmbeddedCanvas(0,fCompositeFrame7353vt,1179,732);
	Int_t wfRootEmbeddedCanvas7363vt = fRootEmbeddedCanvas7363vt->GetCanvasWindowId();
	c2vt = new TCanvas("c2vt", 10, 10, wfRootEmbeddedCanvas7363vt);
	c2vt->Divide(4,4);
	fRootEmbeddedCanvas7363vt->AdoptCanvas(c2vt);
	fCompositeFrame7353vt->AddFrame(fRootEmbeddedCanvas7363vt, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

	// container of "One channel"
	TGCompositeFrame *fCompositeFrame7353vs;
	fCompositeFrame7353vs = fTab683->AddTab("Single CH");
	fCompositeFrame7353vs->SetLayoutManager(new TGVerticalLayout(fCompositeFrame7353vs));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas7363vs = new TRootEmbeddedCanvas(0,fCompositeFrame7353vs,1179,732);
	Int_t wfRootEmbeddedCanvas7363vs = fRootEmbeddedCanvas7363vs->GetCanvasWindowId();
	c2vs = new TCanvas("c2vs", 10, 10, wfRootEmbeddedCanvas7363vs);
	c2vs->Divide(1,1);
	fRootEmbeddedCanvas7363vs->AdoptCanvas(c2vs);
	fCompositeFrame7353vs->AddFrame(fRootEmbeddedCanvas7363vs, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

	// container of "One channel"
	TGCompositeFrame *fCompositeFrame7353;
	fCompositeFrame7353 = fTab683->AddTab("2Fibers@15mm_ScinBar");
	fCompositeFrame7353->SetLayoutManager(new TGVerticalLayout(fCompositeFrame7353));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas7363 = new TRootEmbeddedCanvas(0,fCompositeFrame7353,1179,732);
	Int_t wfRootEmbeddedCanvas7363 = fRootEmbeddedCanvas7363->GetCanvasWindowId();
	c2 = new TCanvas("c2", 10, 10, wfRootEmbeddedCanvas7363);
	c2->Divide(4,2);
	fRootEmbeddedCanvas7363->AdoptCanvas(c2);
	fCompositeFrame7353->AddFrame(fRootEmbeddedCanvas7363, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

	// container of "One channel"
	TGCompositeFrame *fCompositeFrame7353b;
	fCompositeFrame7353b = fTab683->AddTab("2Fibers@10mm_ScinBar");
	fCompositeFrame7353b->SetLayoutManager(new TGVerticalLayout(fCompositeFrame7353b));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas7363b = new TRootEmbeddedCanvas(0,fCompositeFrame7353b,1179,732);
	Int_t wfRootEmbeddedCanvas7363b = fRootEmbeddedCanvas7363b->GetCanvasWindowId();
	c2b = new TCanvas("c2b", 10, 10, wfRootEmbeddedCanvas7363b);
	c2b->Divide(4,2);
	fRootEmbeddedCanvas7363b->AdoptCanvas(c2b);
	fCompositeFrame7353b->AddFrame(fRootEmbeddedCanvas7363b, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

	// container of "TChannel profile"
	TGCompositeFrame *fCompositeFrame7352tt;
	fCompositeFrame7352tt = fTab683->AddTab("Timing");
	fCompositeFrame7352tt->SetLayoutManager(new TGVerticalLayout(fCompositeFrame7352tt));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas7362tt = new TRootEmbeddedCanvas(0,fCompositeFrame7352tt,1179,732);
	Int_t wfRootEmbeddedCanvas7362tt = fRootEmbeddedCanvas7362tt->GetCanvasWindowId();
	myCant = new TCanvas("myCant", 10, 10, wfRootEmbeddedCanvas7362tt);
	myCant->Divide(2,2);
	fRootEmbeddedCanvas7362tt->AdoptCanvas(myCant);
	fCompositeFrame7352tt->AddFrame(fRootEmbeddedCanvas7362tt, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));
//	c2->SetGridx(1);
//	c2->SetGridy(1);

	// container of "TChannel profile"
	TGCompositeFrame *fCompositeFrame7352tts;
	fCompositeFrame7352tts = fTab683->AddTab("FEB Setup");
	fCompositeFrame7352tts->SetLayoutManager(new TGVerticalLayout(fCompositeFrame7352tts));
	// embedded canvas
	TRootEmbeddedCanvas *fRootEmbeddedCanvas7362tts = new TRootEmbeddedCanvas(0,fCompositeFrame7352tts,1179,732);
	Int_t wfRootEmbeddedCanvas7362tts = fRootEmbeddedCanvas7362tts->GetCanvasWindowId();
	myCants = new TCanvas("myCants", 10, 10, wfRootEmbeddedCanvas7362tts);
	myCants->Divide(2,2);
	fRootEmbeddedCanvas7362tts->AdoptCanvas(myCants);
	fCompositeFrame7352tts->AddFrame(fRootEmbeddedCanvas7362tts, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));
//	c2->SetGridx(1);
//	c2->SetGridy(1);


	fTab683->SetTab(0);

	fTab683->Resize(fTab683->GetDefaultSize());
	fHorizontalFrame768->AddFrame(fTab683, new TGLayoutHints(kLHintsRight | kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

	fVerticalFrame734->AddFrame(fHorizontalFrame768, new TGLayoutHints(kLHintsTop | kLHintsExpandX | kLHintsExpandY,2,2,2,2));

	fMainFrame1241->AddFrame(fVerticalFrame734, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY,1,1,1,1));

	fMainFrame1560->AddFrame(fMainFrame1241, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
	fMainFrame1241->MoveResize(0,0,1329,789);

	fMainFrame1314->AddFrame(fMainFrame1560, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
	fMainFrame1560->MoveResize(0,0,1329,789);

	fMainFrame1314->SetMWMHints(kMWMDecorAll,
                        kMWMFuncAll,
                        kMWMInputModeless);
	fMainFrame1314->MapSubwindows();

	fMainFrame1314->Resize(fMainFrame1314->GetDefaultSize());
	fMainFrame1314->MapWindow();
	fMainFrame1314->Resize(1329,789);
	
	
}

TGeoManager *geom;
TEveManager* gEve; 
TEveBoxSet* hits; 
TGFileBrowser *gbrowser;
TFolder* froi;

void DisplayEvent(Bool_t register=kTRUE)
{
	std::cout << " *** DsTau::BuildGeometry() *** " << std::endl;
	gSystem->Load("libGeom");
	
	 if (geom) delete geom;
	// TGeo Geometry
	 geom = new TGeoManager("dstau","DsTau event display");
	Int_t i;
	Float_t z;
	
	//--- define some materials
	TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
	TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);

	//--- define some media
	TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
	TGeoMedium *Al = new TGeoMedium("Aluminium",2, matAl);

	Double_t BoxX = 185.2;
	Double_t BoxY = 185.2;
	Double_t BoxZ = 10;

	Double_t BottomX = 23.;
	Double_t BottomY = 184;
	Double_t BottomZ = 1.5;
	Double_t TopX = 23.;
	Double_t TopY = 184;
	Double_t TopZ = 1;
	
	//--- make the top container volume
	TGeoVolume *top = geom->MakeBox("TOP", Vacuum, BoxX , BoxY, BoxZ);
	geom->SetTopVolume(top);
	
	// Make the elementary assembly of the whole structure
	TGeoVolume *dstau = new TGeoVolumeAssembly("CRTModule");
	
	TGeoBBox *Scin10mm = new TGeoBBox("Scin10mm", BottomX/2 , BottomY/2, BottomZ/2);
	TGeoVolume *vScin10mm = new TGeoVolume("Aluminium",Scin10mm,Al);
	vScin10mm->SetLineColor(kMagenta);
	vScin10mm->SetTransparency(30);
	vScin10mm->SetVisibility(kTRUE);
	
	
	TGeoBBox *Scin15mm = new TGeoBBox("Scin15mm", TopX/2 , TopY/2, TopZ/2);
	TGeoVolume *vScin15mm = new TGeoVolume("Aluminium",Scin15mm,Al);
	vScin15mm->SetLineColor(kMagenta);
	vScin15mm->SetTransparency(30);
	vScin15mm->SetVisibility(kTRUE);
	
	vScin10mm->SetTransparency(30);
	vScin15mm->SetTransparency(30);
	
	Int_t NPlates = 8;
	Double_t CellWidth = 23;
	for( int n = 0; n <NPlates; n++ )
	{
			dstau->AddNode(vScin10mm, n, new TGeoTranslation(0,n*CellWidth,0));
			dstau->AddNode(vScin15mm, n, new TGeoTranslation(0,n*CellWidth,2));
	}
	top->AddNode(dstau, 0, new TGeoTranslation(0,0,0));
	
	//--- close the geometry
	geom->CloseGeometry();

	gStyle->SetPalette(-1);
	TEveManager* evem;
	if (gEve){
		evem = gEve;
		gEve->GetViewers()->DeleteAnnotations();
		gEve->GetCurrentEvent()->DestroyElements();
	} else {
		evem = TEveManager::Create();
	}
	TEveElementList *geometry = new TEveElementList("Geometry");

	// Create Geometry in Event Display
	TGeoNode* node = gGeoManager->GetTopNode();
	TEveGeoTopNode* eveNode = new TEveGeoTopNode(gGeoManager, node);
	eveNode->SetVisLevel(4);
	eveNode->GetNode()->GetVolume()->SetVisibility(kFALSE);
	
	geometry->AddElement(eveNode);
	
	
	gEve->AddGlobalElement(eveNode);

	eveNode->ExpandIntoListTreesRecursively();
//	eveNode->SaveExtract("dstau.root", "Dstau", kFALSE);


	TEveRGBAPalette* pal = new TEveRGBAPalette(0, 100);
	TEveViewer *ev;
	TGLViewer *gv;
	ev = gEve->GetDefaultViewer();
	gEve->GetWindowManager()->HideAllEveDecorations();
	gv = ev->GetGLViewer();
	gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
	gv->SetClearColor(kWhite);

	TEveWindowSlot *slot = 0;
	TEveWindowPack *pack = 0;
	TEveViewer            *T3DView;   
	TEveViewer            *TRPhiView; 
	TEveViewer            *TRhoZView; 
	TEveViewer            *MuonView; 
	slot = TEveWindow::CreateWindowInTab(gEve->GetBrowser()->GetTabRight());
	pack = slot->MakePack(); // slot destroyed
	pack->SetElementName("Multi View");
	pack->SetHorizontal();
	pack->SetShowTitleBar(kFALSE);
	pack->NewSlot()->MakeCurrent(); // new slot created from pack
	T3DView = gEve->SpawnNewViewer("Y-Z View", "");
	T3DView->AddScene(gEve->GetGlobalScene());
	T3DView->AddScene(gEve->GetEventScene());
	gv = T3DView->GetGLViewer();
	gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
	gv->SetClearColor(kWhite);
	pack = pack->NewSlot()->MakePack();
	pack->SetShowTitleBar(kFALSE);
	pack->NewSlot()->MakeCurrent();
	TRPhiView = gEve->SpawnNewViewer("X-Y View", "");
	TRPhiView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
	TRPhiView->AddScene(gEve->GetGlobalScene());
	TRPhiView->AddScene(gEve->GetEventScene());
	gv = TRPhiView->GetGLViewer();
	gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
	gv->SetClearColor(kWhite);
    pack->NewSlot()->MakeCurrent();
	TRhoZView = gEve->SpawnNewViewer("X-Z View", "");
	TRhoZView->GetGLViewer()->SetCurrentCamera(TGLViewer::kCameraOrthoXOZ);
	TRhoZView->AddScene(gEve->GetGlobalScene());
	TRhoZView->AddScene(gEve->GetEventScene());
	gv = TRhoZView->GetGLViewer();
	gv->SetGuideState(TGLUtil::kAxesOrigin, kTRUE, kFALSE, 0);
	gv->SetClearColor(kWhite);

	gEve->GetViewers()->SwitchColorSet();
	gEve->GetDefaultGLViewer()->SetStyle(TGLRnrCtx::kOutline);
	gEve->GetBrowser()->GetTabRight()->SetTab(1);

	gEve->GetBrowser()->StartEmbedding(TRootBrowser::kRight);
	//TRootEmbeddedCanvas  *fEc = new TRootEmbeddedCanvas("ec",);
	gROOT->ProcessLineFast("new TCanvas");
	TCanvas *gcanvas = 0;
	gcanvas = new TCanvas;
	gEve->GetBrowser()->StopEmbedding("Canvas");

	gEve->Redraw3D(kTRUE);

/*
	TString trkname = TString::Format("Track_%d",entry);
	TEveBoxSet* seg = new TEveBoxSet(trkname); 
	seg->Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);
	seg->SetPalette(pal);
	seg->UseSingleColor();
	seg->SetMainColor(kRed);
	seg->SetMainTransparency(1);
	cout << "I am going to extract track segments" << endl;
	//
	// extract segments information from ResultsTTree
	//
	Hits *trk = filereader.tmHits();
	trk->GetHits(entry);
	int charge;
	cout << "The event number for the entry number of " << entry << " is: " << endl;
	cout <<  "    ==> Event: " << trk->Event() << " with the number of hits: " << trk->Nhits() << endl;
	for(Int_t i=0; i<trk->Nhits(); i++){
		double offset = 0.01; // 5 --> 1
		double xminT = trk->PosX(i)/10 - offset;
		double yminT = trk->PosY(i)/10 - offset;
		double zminT = trk->PosZ(i)/10 - offset;
		double size = 0.01;
		//cout << trk->PosX(i)/10 << " " << trk->PosY(i)/10 << " " << trk->PosZ(i)/10 << endl;
		charge = GetCharged(trk->GetPIdCode(i));
		if( charge && trk->PosZ(i)/10 < 10){ // drawing charged tracks only
			seg->AddBox(xminT, yminT, zminT, size, size, size);
			//seg->SetMainColor(kRed);
			//seg->DigitValue(15);
			seg->DigitValue(trk->GetPIdCode(i));
		}
	}

	seg->RefitPlex();
	seg->SetPickable(kTRUE);
	seg->RefMainTrans().SetPos(0,0,0);
	gEve->AddElement(seg);
	gEve->GetBrowser()->GetTabRight()->SetTab(0);
	gEve->Redraw3D(kTRUE);
*/
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Creating TEveBoxSet
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
