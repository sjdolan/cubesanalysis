#include "TreeManager.h"
#include "mppc.h"


#include <iostream>
#include <stdlib.h>
#include <list>
//#include <array>

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif
#ifdef __MAKECINT__
#pragma link off all class;
#pragma link C++ class Mppc+;
#pragma link off all function;
#pragma link off all global;
#pragma link off all typedef;
#endif

using namespace std;

class Mppc;
class TObject;

TreeManager::TreeManager(std::string infilename){
	m_inFilename = infilename;
	cout << "Processing file: " << m_inFilename.c_str() << endl;
	init();
}

TreeManager::~TreeManager(){
	delete  m_treeRes;
	
	delete  m_inFile;
}



void TreeManager::nullify(){
	m_treeRes  = NULL;
	Mppc::releaseThis();
}

//CosmicDisplay *TreeManager::tmCD() {cout << "..tmCD().." << endl; return m_tres;};

void TreeManager::init(){
	nullify();
	m_inFile = new TFile(m_inFilename.c_str(), "READ");
	if( m_inFile->IsOpen() == kFALSE ) return;
	cout << "File opened!" << endl;
	
	m_treeRes = (TTree*) m_inFile->FindObjectAny("mppc");
	if( m_treeRes != NULL ) {
		cout << "m_treeRes found!" << endl;
		m_tres = Mppc::giveThis(m_treeRes,"read");
		if( m_tres != NULL ) {
			cerr << " [1] it is ok!!" << endl;
		}
	}
	
}

