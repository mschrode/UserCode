// $Id: firstDataAnalysis.C,v 1.2 2009/12/29 12:49:21 mschrode Exp $
//
// Plot simple distributions of different samples
// (data, MC,...) from Kalibri ntuples. Meant for
// first data analysis.


#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TVector2.h"


// === Global variables ===
const int maxNJet_ = 50;

TString treeName_ = "DiJetTree";
int nSamples_ = 0;
std::vector<TChain*> chain_;
std::vector<TString> drawOption_;
std::vector<int> markerStyle_;
std::vector<int> color_;

// Distributions after all cuts
std::vector<TH1D*> hPt_;
std::vector<TH1D*> hPtAsym_;
std::vector<TH1D*> hPtCorr_;
std::vector<TH1D*> hEta_;
std::vector<TH1D*> hPhi_;

// N-1 distributions
std::vector<TH1D*> hCutNJets_;
std::vector<TH1D*> hCutVtxNTracks_;
std::vector<TH1D*> hCutVtxPosZ_;
std::vector<TH1D*> hCutRelMet_;
std::vector<TH1D*> hCutJetPt_;
std::vector<TH1D*> hCutJetEta_;
std::vector<TH1D*> hCutDeltaPhi_;
std::vector<TH1D*> hCutEMF_;
std::vector<TH1D*> hCutN90Hits_;
std::vector<TH1D*> hCutFHPD_;
std::vector<TH1D*> hCutFRBX_;

// Cut values
int cutMinNJets_;
int cutMinVtxNTracks_;
double cutMaxVtxPosZ_;
double cutMaxRelMet_;
double cutMinJetPt_;
double cutMaxJetEta_;
double cutMinDeltaPhi_;
double cutMinEMF_;
int cutMinN90Hits_;
double cutMaxFHPD_;
double cutMaxFRBX_;

// Number of events rejected by cuts
std::vector<int> nEvents_;
std::vector<int> nGoodEvents_;
std::vector<int> nCutNJets_;
std::vector<int> nCutLumiBlock_;
std::vector<int> nCutVtx_;
std::vector<int> nCutMet_;
std::vector<int> nCutJetPt_;
std::vector<int> nCutJetEta_;
std::vector<int> nCutJetDeltaPhi_;
std::vector<int> nCutMinEMF_;
std::vector<int> nCutMinN90Hits_;
std::vector<int> nCutMaxFRBX_;
std::vector<int> nCutMaxFHPD_;


bool isInit_ = false;



// === Function declarations ===
void addFile(const TString &fileName, int sample);
void draw(const std::vector<int>& drawnSamples);
void fillHistos();
void init(const TString &treeName, int nSamples);
bool isGoodLumiBlock(int lumiBlockNumber, int runNumber);



// === Main functions ===
void runFirstDataAnalysis() {
  init("DiJetTree",1);
  addFile("/scratch/hh/current/cms/user/mschrode/data/MinBias-BeamCommissioning09-Dec14thReReco_v1/MinBias-BeamCommissioning09-900GeV-Dec14thReReco_v1.root",0);
  fillHistos();

  drawOption_[0] = "PE1";

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
    std::cout << "Plotting distributions... " << std::flush;

    TCanvas *canCutNJets = new TCanvas("canCutNJets","N-1: N jets",0,0,500,500);
    TCanvas *canCutVtxNTracks = new TCanvas("canCutVtxNTracks","N-1: Vtx n tracks",10,20,500,500);
    TCanvas *canCutVtxPosZ = new TCanvas("canCutVtxPosZ","N-1: Vtx z pos",20,40,500,500);
    TCanvas *canCutRelMet = new TCanvas("canCutRelMet","N-1: Rel MET",30,60,500,500);
    TCanvas *canCutJetPt = new TCanvas("canCutJetPt","N-1: Pt",40,80,500,500);
    TCanvas *canCutJetEta = new TCanvas("canCutJetEta","N-1: Eta",50,100,500,500);
    TCanvas *canCutDeltaPhi = new TCanvas("canCutDeltaPhi","N-1: Delta phi",60,120,500,500);
    TCanvas *canCutEMF = new TCanvas("canCutEMF","N-1: EMF",70,140,500,500);
    TCanvas *canCutN90Hits = new TCanvas("canCutN90Hits","N-1: N90 hits",80,160,500,500);
    TCanvas *canCutFHPD = new TCanvas("canCutFHPD","N-1: FHPD",90,180,500,500);
    TCanvas *canCutFRBX = new TCanvas("canCutFRBX","N-1: FRBX",100,200,500,500);

    TCanvas *canPt = new TCanvas("canPt","Pt",0,300,500,500);
    TCanvas *canPtAsym = new TCanvas("canPtAsym","Pt Asymmetry",10,320,500,500);
    TCanvas *canPtCorr = new TCanvas("canPtCorr","Corrected Pt",20,340,500,500);
    TCanvas *canEta = new TCanvas("canEta","Eta",30,360,500,500);
    TCanvas *canPhi = new TCanvas("canPhi","Phi",40,380,500,500);

    std::vector<int>::const_iterator sampleIt = drawnSamples.begin();
    for(; sampleIt != drawnSamples.end(); sampleIt++) {
      TString opt = drawOption_.at(*sampleIt);
      if( sampleIt - drawnSamples.begin() != 0 ) opt += "same";

      canPt->cd();
      hPt_.at(*sampleIt)->Draw(opt);

      canPtAsym->cd();
      hPt_.at(*sampleIt)->Draw(opt);

      canPtCorr->cd();
      hPtCorr_.at(*sampleIt)->Draw(opt);

      canEta->cd();
      hEta_.at(*sampleIt)->Draw(opt);

      canPhi->cd();
      hPhi_.at(*sampleIt)->Draw(opt);
      
      canCutNJets->cd();
      hCutNJets_.at(*sampleIt)->Draw(opt);

      canCutVtxNTracks->cd();
      hCutVtxNTracks_.at(*sampleIt)->Draw(opt);

      canCutVtxPosZ->cd();
      hCutVtxPosZ_.at(*sampleIt)->Draw(opt);

      canCutRelMet->cd();
      hCutRelMet_.at(*sampleIt)->Draw(opt);

      canCutJetPt->cd();
      hCutJetPt_.at(*sampleIt)->Draw(opt);

      canCutJetEta->cd();
      hCutJetEta_.at(*sampleIt)->Draw(opt);
  
      canCutDeltaPhi->cd();
      hCutDeltaPhi_.at(*sampleIt)->Draw(opt);

      canCutEMF->cd();
      hCutEMF_.at(*sampleIt)->Draw(opt);

      canCutN90Hits->cd();
      hCutN90Hits_.at(*sampleIt)->Draw(opt);

      canCutFHPD->cd();
      hCutFHPD_.at(*sampleIt)->Draw(opt);
	
      canCutFRBX->cd();
      hCutFRBX_.at(*sampleIt)->Draw(opt);
    }
    std::cout << "ok\n";
  } else {
    std::cerr << "ERROR: Objects not initialized. Run init() first.\n";
  }
}


void fillHistos() {
  if( isInit_ ) {
    std::cout << "Selecting events and filling distributions... " << std::flush;
    
    // Init read quantities
    unsigned int runNumber = 0;
    unsigned int lumiBlockNumber = 0;
    int vtxNTracks = 0;
    float vtxPosZ = 0.;
    float sumEt = 0.;
    float met = 0.;
    int nObjJet = 0;
    float jetPt[maxNJet_];
    float jetEta[maxNJet_];
    float jetPhi[maxNJet_];
    float jetEMF[maxNJet_];
    int jetN90Hits[maxNJet_];
    float jetFHPD[maxNJet_];
    float jetFRBX[maxNJet_];
    float jetCorrL2L3[maxNJet_];

    // Loop over samples
    for(int i = 0; i < nSamples_; i++) {

      // Reset histos
      hPt_[i]->Reset();

      // Set branch addresses
      chain_[i]->SetBranchAddress("RunNumber",&runNumber);
      chain_[i]->SetBranchAddress("LumiBlockNumber",&lumiBlockNumber);
      chain_[i]->SetBranchAddress("VtxNTracks",&vtxNTracks);
      chain_[i]->SetBranchAddress("VtxPosZ",&vtxPosZ);
      chain_[i]->SetBranchAddress("Met",&met);
      chain_[i]->SetBranchAddress("MetSum",&sumEt);
      chain_[i]->SetBranchAddress("NobjJet",&nObjJet);
      chain_[i]->SetBranchAddress("JetPt",jetPt);
      chain_[i]->SetBranchAddress("JetEta",jetEta);
      chain_[i]->SetBranchAddress("JetPhi",jetPhi);
      chain_[i]->SetBranchAddress("JetEMF",jetEMF);
      chain_[i]->SetBranchAddress("JetN90Hits",jetN90Hits);
      chain_[i]->SetBranchAddress("JetFHPD",jetFHPD);
      chain_[i]->SetBranchAddress("JetFRBX",jetFRBX);
      chain_[i]->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  
      // Loop over tree entries
      for(int n = 0; n < 100; n++ ) {//chain_[i]->GetEntries(); n++) {
	chain_[i]->GetEntry(n);

	if( nObjJet > maxNJet_ ) {
	  std::cerr << "WARNING: nObjJet = " << nObjJet << " > maxNJet_. Skipping event.\n";
	  continue;
	}

	// Event selection
	nEvents_[i]++;
	bool isGood = true;

	// Cut on number of jets
	hCutNJets_[i]->Fill(nObjJet);
	if( nObjJet < cutMinNJets_ ) {
	  nCutNJets_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on luminosity block
	if( !isGoodLumiBlock(lumiBlockNumber,runNumber) ) {
	  nCutLumiBlock_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on primary vertex
	hCutVtxNTracks_[i]->Fill(vtxNTracks);
	hCutVtxPosZ_[i]->Fill(vtxPosZ);
	if( vtxNTracks < cutMinVtxNTracks_ || std::abs(vtxPosZ) > cutMaxVtxPosZ_ ) {
	  nCutVtx_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on MET
	double relMet = met/sumEt;
	hCutRelMet_[i]->Fill(relMet);
	if( relMet > cutMaxRelMet_ ) {
	  nCutMet_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet pt
	hCutJetPt_[i]->Fill(jetPt[1]);
	if( jetPt[1] < cutMinJetPt_ ) {
	  nCutJetPt_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet eta
	hCutJetEta_[i]->Fill(jetEta[0]);
	hCutJetEta_[i]->Fill(jetEta[1]);
	if( std::abs(jetEta[0]) > cutMaxJetEta_ || std::abs(jetEta[1]) > cutMaxJetEta_ ) {
	  nCutJetEta_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet Delta phi
	double deltaPhi = TVector2::Phi_mpi_pi(jetPhi[0] - jetPhi[1]);
	hCutDeltaPhi_[i]->Fill(deltaPhi);
	if( std::abs(deltaPhi) < cutMinDeltaPhi_ ) {
	  nCutJetDeltaPhi_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet EMF
	hCutEMF_[i]->Fill(jetEMF[0]);
	hCutEMF_[i]->Fill(jetEMF[1]);
	if( jetEMF[0] < cutMinEMF_ || jetEMF[1] < cutMinEMF_) {
	  nCutMinEMF_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet n90Hits
	hCutN90Hits_[i]->Fill(jetN90Hits[0]);
	hCutN90Hits_[i]->Fill(jetN90Hits[1]);
	if( jetN90Hits[0] < cutMinN90Hits_ || jetN90Hits[1] < cutMinN90Hits_ ) {
	  nCutMinN90Hits_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet FRBX
	hCutFRBX_[i]->Fill(jetFRBX[0]);
	hCutFRBX_[i]->Fill(jetFRBX[1]);
	if( jetFRBX[0] > cutMaxFRBX_ || jetFRBX[1] > cutMaxFRBX_ ) {
	  nCutMaxFRBX_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet FHPD
	hCutFHPD_[i]->Fill(jetFHPD[0]);
	hCutFHPD_[i]->Fill(jetFHPD[1]);
	if( jetFHPD[0] > cutMaxFHPD_ || jetFHPD[1] > cutMaxFHPD_ ) {
	  nCutMaxFHPD_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;


	// Fill distributions of selected events
	nGoodEvents_[i]++;
	hPtAsym_[i]->Fill((jetPt[0]-jetPt[1])/(jetPt[0]+jetPt[1]));
	for(int j = 0; j < 2; j++) {
	  hPt_[i]->Fill(jetPt[j]);
	  hPtCorr_[i]->Fill(jetCorrL2L3[j]*jetPt[j]);
	  hEta_[i]->Fill(jetEta[j]);
	  hPhi_[i]->Fill(jetPhi[j]);
	}
      } // End of loop over tree entries
    } // End of loop over samples
    std::cout << "ok\n";
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

    drawOption_ = std::vector<TString>(nSamples_);
    markerStyle_ = std::vector<int>(nSamples_);
    color_ = std::vector<int>(nSamples_);

    // Cuts
    cutMinVtxNTracks_ = 1;
    cutMaxVtxPosZ_ = 20.;
    cutMinNJets_ = 2;
    cutMaxRelMet_ = 0.5;
    cutMinJetPt_ = 10.;
    cutMaxJetEta_ = 3.;
    cutMinDeltaPhi_ = 2.1;
    cutMinEMF_ = 0.01;
    cutMinN90Hits_ = 2;
    cutMaxFHPD_ = 0.98;
    cutMaxFRBX_ = 0.98;

    nEvents_ = std::vector<int>(nSamples_);
    nGoodEvents_ = std::vector<int>(nSamples_);
    nCutNJets_ = std::vector<int>(nSamples_);
    nCutLumiBlock_ = std::vector<int>(nSamples_);
    nCutVtx_ = std::vector<int>(nSamples_);
    nCutMet_ = std::vector<int>(nSamples_);
    nCutJetPt_ = std::vector<int>(nSamples_);
    nCutJetEta_ = std::vector<int>(nSamples_);
    nCutJetDeltaPhi_ = std::vector<int>(nSamples_);
    nCutMinEMF_ = std::vector<int>(nSamples_);
    nCutMinN90Hits_ = std::vector<int>(nSamples_);
    nCutMaxFRBX_ = std::vector<int>(nSamples_);
    nCutMaxFHPD_ = std::vector<int>(nSamples_);
    for(int i = 0; i < nSamples_; i++) {
      nEvents_[i] = 0;
      nGoodEvents_[i] = 0;
      nCutNJets_[i] = 0;
      nCutLumiBlock_[i] = 0;
      nCutVtx_[i] = 0;
      nCutMet_[i] = 0;
      nCutJetPt_[i] = 0;
      nCutJetEta_[i] = 0;
      nCutJetDeltaPhi_[i] = 0;
      nCutMinEMF_[i] = 0;
      nCutMinN90Hits_[i] = 0;
      nCutMaxFRBX_[i] = 0;
      nCutMaxFHPD_[i] = 0;
    }

    // Histos
    chain_ = std::vector<TChain*>(nSamples_);
    hPt_ = std::vector<TH1D*>(nSamples_);
    hPtAsym_ = std::vector<TH1D*>(nSamples_);
    hPtCorr_ = std::vector<TH1D*>(nSamples_);
    hEta_ = std::vector<TH1D*>(nSamples_);
    hPhi_ = std::vector<TH1D*>(nSamples_);

    // N-1 distributions
    hCutNJets_ = std::vector<TH1D*>(nSamples_);
    hCutVtxNTracks_ = std::vector<TH1D*>(nSamples_);
    hCutVtxPosZ_ = std::vector<TH1D*>(nSamples_);
    hCutRelMet_ = std::vector<TH1D*>(nSamples_);
    hCutJetPt_ = std::vector<TH1D*>(nSamples_);
    hCutJetEta_ = std::vector<TH1D*>(nSamples_);
    hCutDeltaPhi_ = std::vector<TH1D*>(nSamples_);
    hCutEMF_ = std::vector<TH1D*>(nSamples_);
    hCutN90Hits_ = std::vector<TH1D*>(nSamples_);
    hCutFHPD_ = std::vector<TH1D*>(nSamples_);
    hCutFRBX_ = std::vector<TH1D*>(nSamples_);

    TString name;
    for(int i = 0; i < nSamples_; i++) {
      chain_[i] = new TChain(treeName_);
      
      drawOption_[i] = "HIST";
      markerStyle_[i] = 20;
      color_[i] = 1;
      
      name = "hPt";
      name += i;
      hPt_[i] = new TH1D(name,"Selected dijet events",40,0,100);
      hPt_[i]->GetXaxis()->SetTitle("p_{T} (GeV)");
      hPt_[i]->GetYaxis()->SetTitle("Number of jets");
      hPt_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hPtAsym";
      name += i;
      hPtAsym_[i] = new TH1D(name,"Selected dijet events",40,-1.5,1.5);
      hPtAsym_[i]->GetXaxis()->SetTitle("p_{T} asymmetry");
      hPtAsym_[i]->GetYaxis()->SetTitle("Number of dijet events");
      hPtAsym_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hPtCorr";
      name += i;
      hPtCorr_[i] = new TH1D(name,"Selected dijet events",40,-1.5,1.5);
      hPtCorr_[i]->GetXaxis()->SetTitle("L2L3 corrected p_{T} (GeV)");
      hPtCorr_[i]->GetYaxis()->SetTitle("Number of jets");
      hPtCorr_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hEta";
      name += i;
      hEta_[i] = new TH1D(name,"Selected dijet events",40,-5.,5.);
      hEta_[i]->GetXaxis()->SetTitle("#eta");
      hEta_[i]->GetYaxis()->SetTitle("Number of jets");
      hEta_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hPhi";
      name += i;
      hPhi_[i] = new TH1D(name,"Selected dijet events",40,-3.2,3.2);
      hPhi_[i]->GetXaxis()->SetTitle("#phi");
      hPhi_[i]->GetYaxis()->SetTitle("Number of jets");
      hPhi_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutNJets";
      hCutNJets_[i] = new TH1D(name,"n - 1 plot",30,0,30);
      hCutNJets_[i]->GetXaxis()->SetTitle("Number of jets");
      hCutNJets_[i]->GetYaxis()->SetTitle("Number of events");
      hCutNJets_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutVtxNTracks";
      hCutVtxNTracks_[i] = new TH1D(name,"n - 1 plot",20,0,20);
      hCutVtxNTracks_[i]->GetXaxis()->SetTitle("Number of primary vertex tracks");
      hCutVtxNTracks_[i]->GetYaxis()->SetTitle("Number of events");
      hCutVtxNTracks_[i]->SetMarkerStyle(markerStyle_[i]);
       
      name = "hCutVtxPosZ";
      hCutVtxPosZ_[i] = new TH1D(name,"n - 1 plot",40,-50,50);
      hCutVtxPosZ_[i]->GetXaxis()->SetTitle("Primary vertex z-position (cm)");
      hCutVtxPosZ_[i]->GetYaxis()->SetTitle("Number of events");
      hCutVtxPosZ_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutRelMet";
      hCutRelMet_[i] = new TH1D(name,"n - 1 plot",40,0,1.1);
      hCutRelMet_[i]->GetXaxis()->SetTitle("MET / sumET");
      hCutRelMet_[i]->GetYaxis()->SetTitle("Number of events");
      hCutRelMet_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutJetPt";
      hCutJetPt_[i] = new TH1D(name,"n - 1 plot",40,0,100);
      hCutJetPt_[i]->GetXaxis()->SetTitle("p^{jet2}_{T} (GeV)");
      hCutJetPt_[i]->GetYaxis()->SetTitle("Number of jets");
      hCutJetPt_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutJetEta";
      hCutJetEta_[i] = new TH1D(name,"n - 1 plot",40,-5.2,5.2);
      hCutJetEta_[i]->GetXaxis()->SetTitle("#eta");
      hCutJetEta_[i]->GetYaxis()->SetTitle("Number of jets");
      hCutJetEta_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutDeltaPhi";
      hCutDeltaPhi_[i] = new TH1D(name,"n - 1 plot",40,-3.2,3.2);
      hCutDeltaPhi_[i]->GetXaxis()->SetTitle("#Delta#phi(jet1,jet2)");
      hCutDeltaPhi_[i]->GetYaxis()->SetTitle("Number of events");
      hCutDeltaPhi_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutEMF";
      hCutEMF_[i] = new TH1D(name,"n - 1 plot",40,-2,1);
      hCutEMF_[i]->GetXaxis()->SetTitle("f_{EM}");
      hCutEMF_[i]->GetYaxis()->SetTitle("Number of jets");
      hCutEMF_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutN90Hits";
      hCutN90Hits_[i] = new TH1D(name,"n - 1 plot",20,0,20);
      hCutN90Hits_[i]->GetXaxis()->SetTitle("n90Hits");
      hCutN90Hits_[i]->GetYaxis()->SetTitle("Number of jets");
      hCutN90Hits_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutFHPD";
      hCutFHPD_[i] = new TH1D(name,"n - 1 plot",40,0,1);
      hCutFHPD_[i]->GetXaxis()->SetTitle("f_{HPD}");
      hCutFHPD_[i]->GetYaxis()->SetTitle("Number of jets");
      hCutFHPD_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hCutFRBX";
      hCutFRBX_[i] = new TH1D(name,"n - 1 plot",40,0,1);
      hCutFRBX_[i]->GetXaxis()->SetTitle("f_{RBX}");
      hCutFRBX_[i]->GetYaxis()->SetTitle("Number of jets");
      hCutFRBX_[i]->SetMarkerStyle(markerStyle_[i]);
    }

    isInit_ = true;

    std::cout << "ok\n";
  }
}


bool isGoodLumiBlock(int lumiBlockNumber, int runNumber) {
  bool isGood = true;
  if( runNumber == 123596 ) {
    if( lumiBlockNumber < 69 ) isGood = false;
  }

  return isGood;
}
