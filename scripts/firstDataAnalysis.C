#include "TCanvas.h"
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

const int maxNJet_ = 50;
const int nSamples_ = 2;

bool isInit_ = false;

// -1: All
// 0: MC
// 1: Data
int plottedSamples_ = 0; 
TString treeName_ = "DiJetTree";

// 0: MC
// 1: Data
TChain *chain_[nSamples_];

TH1F *hPt_[nSamples_];


void init(const TString &treeName, int plottedSamples) {
  if( isInit_ ) {
    std::cout << "Objects already initialized. Skipping init().\n";
  } else {
    plottedSamples_ = plottedSamples;
    treeName_ = treeName;

    TString name;
    for(int i = 0; i < nSamples_; i++) {
      chain_[i] = new TChain(treeName_);
      
      name = "hPt";
      name += i;
      hPt_[i] = new TH1F(name,"",40,0,40);
      hPt_[i]->GetXaxis()->SetTitle("p_{T} (GeV)");
      hPt_[i]->GetXaxis()->SetTitle("Number of jets");
    }

    isInit_ = true;
  }
}


void fillHistos(bool isData) {
  if( isInit_ ) {
    int idx = 0;
    if( isData ) idx = 1;

    // Reset histos
    for(int i = 0; i < nSamples_; i++) {
      hPt_[idx]->Reset();
    }

    int nObjJet = 0;
    float jetPt[maxNJet_];

    chain_[idx]->SetBranchAddress("NobjJet",&nObjJet);
    chain_[idx]->SetBranchAddress("JetPt",jetPt);
  
    for(int n = 0; n < chain_[idx]->GetEntries(); n++) {
      chain_[idx]->GetEntry(n);

      if( nObjJet > maxNJet_ ) continue;

      hPt_[idx]->Fill(jetPt[0]);
    }
  } else {
    std::cerr << "ERROR: Objects not initialized. Run init() first.\n";
  }
}


void draw() {
  if( isInit_ ) {
    TCanvas *canPt = new TCanvas("canPt","Jet Pt",500,500);
    
    bool isFirst = true;
    for(int i = 0; i < nSamples_; i++) {
      if( plottedSamples_ == i || plottedSamples_ == -1 ) {
	canPt->cd();
	if( isFirst ) {
	  hPt_[i]->Draw();
	  isFirst = false;
	} else {
	  hPt_[i]->Draw("same");
	}
      }
    }
  } else {
    std::cerr << "ERROR: Objects not initialized. Run init() first.\n";
  }
}


void addFile(const TString &fileName, bool isData) {
  if( isInit_ ) {
    if( isData ) {
      chain_[0].AddFile(fileName);
    } else {
      chain_[1].AddFile(fileName);
    }
  } else {
    std::cerr << "ERROR: Objects not initialized. Run init() first.\n";
  }
}

