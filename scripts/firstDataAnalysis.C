// $Id: $
//
// Plot simple distributions of different samples
// (data, MC,...) from Kalibri ntuples. Meant for
// first data analysis.


#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"

// === Global variables ===
const int maxNJet_ = 50;

TString treeName_ = "DiJetTree";
int nSamples_ = 0;
std::vector<TChain*> chain_;
std::vector<TH1F*> hPt_;

bool isInit_ = false;



// === Function declarations ===
void addFile(const TString &fileName, int sample);
void draw(const std::vector<int>& drawnSamples);
void fillHistos();
void init(const TString &treeName, int nSamples);



// === Main functions ===
void runFirstDataAnalysis() {
  init("DiJetTree",1);
  addFile("/scratch/hh/current/cms/user/mschrode/data/MinBias-BeamCommissioning09-Dec14thReReco_v1/MinBias-BeamCommissioning09-900GeV-Dec14thReReco_v1.root",0);
  fillHistos();

  std::vector<int> drawnSamples(1);
  drawnSamples.at(0) = 0;
  draw(drawnSamples);
}



// === Function implemenations ===
void addFile(const TString &fileName, int sample) {
  if( isInit_ ) {
    if( sample >=0 && sample < nSamples_ ) {
      std::cout << "Adding file '" << fileName << "' to sample '" << sample << "'... " << std::flush;
      chain_.at(sample)->AddFile(fileName);
      std::cout << "ok\n";
    } else {
      std::cerr << "ERROR: There is no sample with index '" << sample << "'. Skipping.\n";
    }
  } else {
    std::cerr << "ERROR: Objects not initialized. Run init() first.\n";
  }
}


void draw(const std::vector<int>& drawnSamples) {
  if( isInit_ ) {
    TCanvas *canPt = new TCanvas("canPt","Jet Pt",500,500);
    
    std::vector<int>::const_iterator sampleIt = drawnSamples.begin();
    for(; sampleIt != drawnSamples.end(); sampleIt++) {
      if( sampleIt - drawnSamples.begin() == 0 ) {
	canPt->cd();
	hPt_.at(*sampleIt)->Draw();
      } else {
	hPt_.at(*sampleIt)->Draw("same");
      }
    }
  } else {
    std::cerr << "ERROR: Objects not initialized. Run init() first.\n";
  }
}


void fillHistos() {
  if( isInit_ ) {
    // Init read quantities
    int nObjJet = 0;
    float jetPt[maxNJet_];

    // Loop over samples
    for(int i = 0; i < nSamples_; i++) {

      // Reset histos
      hPt_[i]->Reset();

      // Set branch addresses
      chain_[i]->SetBranchAddress("NobjJet",&nObjJet);
      chain_[i]->SetBranchAddress("JetPt",jetPt);
  
      // Loop over tree entries
      for(int n = 0; n < chain_[i]->GetEntries(); n++) {
	chain_[i]->GetEntry(n);

	if( nObjJet > maxNJet_ ) continue;

	hPt_[i]->Fill(jetPt[0]);
      } // End of loop over tree entries
    } // End of loop over samples
  } else {
    std::cerr << "ERROR: Objects not initialized. Run init() first.\n";
  }
}


void init(const TString &treeName, int nSamples) {
  if( isInit_ ) {
    std::cout << "Objects already initialized. Skipping init().\n";
  } else {
    std::cout << "Initializing objects... " << std::flush;

    treeName_ = treeName;

    nSamples_ = nSamples;
    chain_ = std::vector<TChain*>(nSamples_);
    hPt_ = std::vector<TH1F*>(nSamples_);

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

    std::cout << "ok\n";
  }
}



