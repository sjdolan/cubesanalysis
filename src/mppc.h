#ifndef MPPC_H
#define MPPC_H

#define DEBUG_rg 0 

#include <TROOT.h>
#include <TTree.h>
#include <TBranch.h>
#include <TNamed.h>
#include <TObject.h>
#include <TFile.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <TTreePerfStats.h>
//#include <TStopwatch.h>

using namespace std;

class Mppc : public TObject {
 public :
  static const Int_t  channel =32;
  UChar_t         mac5;
  //UChar_t         flags;
  UShort_t        chg[32];
  UInt_t          ts0;
  UInt_t          ts1;
  UInt_t          ts0_ref;
  UInt_t          ts1_ref;
  UInt_t          coincidence;
  
  // List of branches
  TBranch        *b_mac5;   //!
  //TBranch        *b_flags;  //!
  TBranch        *b_chg;   //!
  TBranch        *b_ts0;   //!
  TBranch        *b_ts1;   //!
  TBranch        *b_ts0_ref;   //!
  TBranch        *b_ts1_ref;   //!
  TBranch        *b_coinc;
  //
  Mppc();
		~Mppc(){};
		static Mppc* giveThis();
		static Mppc* giveThis(TTree* aTree, const std::string& option);
		static void releaseThis();
		void setupRead(const std::string& option = "");
		bool SetBranchAddresses();
		Mppc(TTree* aTree, const std::string& option);
		TTree* GetInputTree();
		void GetMppc(Long64_t entry, bool isSFGDChannelSwap);
		
		UChar_t      Mac()                        { return m_mac;        };
		//UChar_t      Flag()                       { return m_flag;       };
		UInt_t       TS0()                        { return m_ts0;        };
		UInt_t       TS0Ref()                     { return m_ts1;        };
		UInt_t       TS1()                        { return m_ts0ref;     };
		UInt_t       TS1Ref()                     { return m_ts1ref;     };
		UInt_t       ADC(Int_t index )            { return m_chg[index]; };
		UInt_t       Coincidence()                { return m_coin;       }
	private:
		TTree* m_treeIn;
		UChar_t m_mac, m_flag;
		UInt_t 	m_ts0, m_ts0ref, m_ts1, m_ts1ref; 
		UInt_t m_chg[channel], m_coin;
		static Mppc* m_instance;
		ClassDef(Mppc,0)
	};
	// initialisation of the Mppc pointer
	Mppc* Mppc::m_instance = 0;
	Mppc::Mppc(){	}
	// to get a unique instance 
	inline Mppc* Mppc::giveThis()
	{
		if(DEBUG_rg) cout << "debug::m_instance Mppc " << m_instance << endl;
		if (0 == m_instance){
			cout << "Mppc::giveThis error not constructed properly " << endl;
		}
		return m_instance;
	}
	// Open a TTree for reading (writing option not now)
	inline Mppc* Mppc::giveThis(TTree* aTree, const std::string& option)
	{
		if(DEBUG_rg) cout << "debug::m_instance Mppc " << m_instance << endl;
		if (0 == m_instance){
			m_instance = new Mppc(aTree, option);
		} else{
			cout << "Mppc::giveThis Warning " << aTree->GetTitle() << endl;
		}
		return m_instance;
	}
	// Delete unique instance
	inline void Mppc::releaseThis() {
		if ( m_instance != 0 ) {
			delete m_instance;
			m_instance = 0;
		}
	}
	// constructor for one TTree with read/write option
	Mppc::Mppc(TTree* aTree, const std::string& option){
		if(option=="read"){
			Long64_t nevent = 9999999;
			m_treeIn = aTree;
			if(DEBUG_rg) cout << " ....before setupread.... " << endl;
//			Int_t cachesize = 10000000;
			setupRead();
//			Long64_t nentries = m_treeIn->GetEntries();
//			nevent = TMath::Min(nevent,nentries);
//			m_treeIn->SetCacheSize(cachesize);
//			m_treeIn->SetCacheLearnEntries(1);
//			m_treeIn->SetCacheEntryRange(0,nevent);
		}
	}
	// setup the input tree for reading
	void Mppc::setupRead(const std::string& option){
		if(DEBUG_rg) cout << " ....inside setupread.... " << endl;
		if(DEBUG_rg) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
		if(!m_treeIn){
			cout << "setupRead error: m_treeIn undefined " << endl;
			exit(1);
		}
		if(SetBranchAddresses()){}else{
			cerr << "TBranch error.."<< endl;
		}
	}
	// Setting up Branch content of input TTree
	bool Mppc::SetBranchAddresses(){
		if(DEBUG_rg) cout << " ....inside setbranchaddresses.... " << endl;
		m_treeIn->SetMakeClass(1);
		m_treeIn->SetBranchAddress("mac5", &mac5, &b_mac5);
		//m_treeIn->SetBranchAddress("flags", &flags, &b_flags);
		m_treeIn->SetBranchAddress("chg", chg, &b_chg);
		m_treeIn->SetBranchAddress("ts0", &ts0, &b_ts0);
		m_treeIn->SetBranchAddress("ts1", &ts1, &b_ts1);
		m_treeIn->SetBranchAddress("ts0_ref", &ts0_ref, &b_ts0_ref);
		m_treeIn->SetBranchAddress("ts1_ref", &ts1_ref, &b_ts1_ref);
		m_treeIn->SetBranchAddress("coincidence", &coincidence, &b_coinc);
		//
		if(DEBUG_rg) cout << " m_treeIn::" << m_treeIn->GetEntries() << endl;
		return true;
	}
	// accessors for the input TTree
	TTree* Mppc::GetInputTree() {return m_treeIn;};
	// get all branch contents of input TTree for entry i
	void Mppc::GetMppc(Long64_t entry, bool isSFGDChannelSwap=true)
	{
		if (!m_treeIn) {
			cout << "Mppc::getEntry error" << endl;
			exit(1);
		}
		//
		m_treeIn->GetEntry(entry);
		if(DEBUG_rg) m_treeIn->Show(entry);
		m_mac    = mac5;
		//m_flag   = flags;
		m_ts0    = ts0;
		m_ts0ref = ts0_ref;
		m_ts1    = ts1;
		m_ts1ref = ts1_ref;
		m_coin   = coincidence;
		for(int i=0; i<32;i++)
		{
			// Correct for channel swap for SFGD run:
			if (i==8 && isSFGDChannelSwap) m_chg[8] = chg[18];
			else if (i==18 && isSFGDChannelSwap) m_chg[18] = chg[8];
			else m_chg[i] = chg[i];
		}
		if(DEBUG_rg) m_treeIn->Print();
	}



#endif

