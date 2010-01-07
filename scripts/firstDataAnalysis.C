// $Id: firstDataAnalysis.C,v 1.4 2009/12/29 19:05:45 mschrode Exp $
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
#include "TLegend.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
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
TString legendHeader_;
TRandom3 *rand_;

int nCuts_;
std::vector< std::vector<TH1D*> > hCuts_;    // N-1 distributions
std::vector<double> cutValue_;               // Cut values
std::vector< std::vector<int> > nCutEvents_; // Number of events rejected by cuts
std::vector<int> nEvents_;
std::vector<int> nGoodEvents_;
std::vector<int> nCutLumiBlock_;


// Distributions after all cuts
std::vector<TH1D*> hPt_;
std::vector<TH1D*> hPtAsym_;
std::vector<TH1D*> hPtCorr_;
std::vector<TH1D*> hEta_;
std::vector<TH1D*> hPhi_;
std::vector<TH1D*> hDeltaPhi_;

bool isInit_ = false;



// === Function declarations ===
void addFile(const TString &fileName, int sample);
void draw(int normSample, const TString &sampleName);
void fillHistos(int events = -1);
void init(const TString &treeName, int nSamples);
bool isGoodLumiBlock(int lumiBlockNumber, int runNumber);
void normaliseDistributions(int normSample);
void printCutFlow(int normSample);



// === Main functions ===
void firstDataAnalysis(int nEvts, const TString &sample = "900GeV") {
  init("DiJetTree",2);
  TString dir = "/scratch/hh/current/cms/user/mschrode/";
  if( sample == "900GeV" ) {
    TString name = dir;
    name += "data/MinBias-BeamCommissioning09-Dec14thReReco_v1/MinBias-BeamCommissioning09-900GeV-Dec14thReReco_v1.root";
    addFile(name,0);
    for(int f = 1; f <= 3; f++) {
      name = dir;
      name += "mc/MinBias-Summer09-DESIGN_3X_V8A_900GeV-v1/MinBias-Summer09-DESIGN_3X_V8A_900GeV-v1__";
      name += f;
      name += ".root";
      addFile(name,1);
    }
    legendHeader_ = "900 GeV MinBias";
  } else if( sample == "2360GeV" ) {
    TString name = dir;
    name += "data/MinBias-BeamCommissioning09-Dec14thReReco_v1/MinBias-BeamCommissioning09-2360GeV-Dec14thReReco_v1.root";
    addFile(name,0);
    for(int f = 1; f <= 5; f++) {
      name = dir;
      name += "mc/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1__";
      name += f;
      name += ".root";
      addFile(name,1);
    }
    legendHeader_ = "2360 GeV MinBias";
  }
  
  fillHistos(nEvts);
  normaliseDistributions(0);

  drawOption_[0] = "PE1";
  draw(0,sample);

  printCutFlow(0);
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


void draw(int normSample, const TString &sampleName) {
  if( isInit_ ) {
    std::cout << "Plotting distributions... " << std::flush;

    gStyle->SetOptStat(0);
    std::vector<TCanvas*> canCuts(nCuts_);
    for(int c = 0; c < nCuts_; c++) {
      TString name = "canCut";
      name += c;
      TString title = "Cut_";
      title += c;
      canCuts.at(c) = new TCanvas(name,title,10*c,20*c,500,500);
    }
    std::vector<TCanvas*> canDist(6);
    canDist[0] = new TCanvas("canPt","Pt",0,300,500,500);
    canDist[1] = new TCanvas("canPtAsym","PtAsymmetry",10,320,500,500);
    canDist[2] = new TCanvas("canPtCorr","CorrectedPt",20,340,500,500);
    canDist[3] = new TCanvas("canEta","Eta",30,360,500,500);
    canDist[4] = new TCanvas("canPhi","Phi",40,380,500,500);
    canDist[5] = new TCanvas("canDeltaPhi","DeltaPhi",50,400,500,500);

    std::vector<TLegend*> leg(nCuts_+6);
    for(size_t i = 0; i < leg.size(); i++) {
      leg[i] = new TLegend(0.5,0.6,0.93,0.85);
      leg[i]->SetBorderSize(0);
      leg[i]->SetFillColor(0);
      leg[i]->SetTextFont(42);
      leg[i]->SetHeader(legendHeader_);
    }
    std::vector<bool> logY(nCuts_+6,false);
    logY[0] = true;
    logY[2] = true;
    logY[6] = true;
    logY[10] = true;

    for(int i = 0; i < nSamples_; i++) {
      TString opt = drawOption_.at(i);
      if( i ) opt += "same";

      for(int n = 0; n < 6; n++) {
	TH1D *h = hPt_[i];
	if( n == 1 ) h = hPtAsym_[i];
	else if( n == 2 ) h = hPtCorr_[i];
	else if( n == 3 ) h = hEta_[i];
	else if( n == 4 ) h = hPhi_[i];
	else if( n == 5 ) h = hDeltaPhi_[i];

	if( i == normSample ) {
	  char text[50];
      	  sprintf(text,"Data: %i entries",static_cast<int>(h->GetEntries()));
	  leg[n]->AddEntry(h,text,"PL");
	} else {
	  leg[n]->AddEntry(h,"MC: normalized","L");
	}

	if( i == 0 ) {
	  double max = h->GetMaximum();
	  h->GetYaxis()->SetRangeUser(0.,2.*max);
	  if( logY[n] ) {
	    h->GetYaxis()->SetRangeUser(0.1,5*h->GetMaximum());
	  }
	}

	canDist[n]->cd();
	h->Draw(opt);

	if( i == nSamples_ - 1 ) {
	  if( logY[n] ) canDist[n]->SetLogy(1);
	  leg[n]->Draw("same");
	  TString name = "plots/";
	  name += sampleName;
	  name += "_";
	  name += canDist[n]->GetTitle();
	  name += ".eps";
	  canDist[n]->SaveAs(name,"eps");
	}
      }      

      for(int c = 0; c < nCuts_; c++) {
	canCuts.at(c)->cd();
	TH1D *h = hCuts_[c][i];
	int idx = c + 6;
	if( i == normSample ) {
	  char text[50];
      	  sprintf(text,"Data: %i entries",static_cast<int>(h->GetEntries()));
	  leg[idx]->AddEntry(h,text,"PL");
	} else {
	  leg[idx]->AddEntry(h,"MC: normalized","L");
	}
	if( i == 0 ) {
	  h->GetYaxis()->SetRangeUser(0.,1.5*h->GetMaximum());
	  if( logY[idx] ) {
	    h->GetYaxis()->SetRangeUser(0.1,5*h->GetMaximum());
	  }
	}
	h->Draw(opt);
	if( i == nSamples_ - 1 ) {
	  if( logY[idx] ) canCuts.at(c)->SetLogy(1); 
	  leg[idx]->Draw("same");
	  TString name = "plots/";
	  name += sampleName;
	  name += "_";
	  name += canCuts[c]->GetTitle();
	  name += ".eps";
	  canCuts[c]->SaveAs(name,"eps");
	}
      }
    }
    std::cout << "ok\n";
  } else {
    std::cerr << "ERROR: Objects not initialized. Run init() first.\n";
  }
}


void fillHistos(int events) {
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
      int nEntries = chain_[i]->GetEntries();
      if( events > 0 && nEntries > events ) nEntries = events;
      for(int n = 0; n < nEntries; n++) {
	chain_[i]->GetEntry(n);

	if( nObjJet > maxNJet_ ) {
	  std::cerr << "WARNING: nObjJet = " << nObjJet << " > maxNJet_. Skipping event.\n";
	  continue;
	}

	// Event selection
	nEvents_[i]++;
	bool isGood = true;

	// Cut on luminosity block
	if( !isGoodLumiBlock(lumiBlockNumber,runNumber) ) {
	  nCutLumiBlock_[i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on number of jets
	hCuts_[0][i]->Fill(nObjJet);
	if( nObjJet < cutValue_[0] ) {
	  nCutEvents_[0][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on primary vertex
	hCuts_[1][i]->Fill(vtxNTracks);
	if( vtxNTracks < cutValue_[1] ) {
	  nCutEvents_[1][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	hCuts_[2][i]->Fill(vtxPosZ);
	if( std::abs(vtxPosZ) > cutValue_[2] ) {
	  nCutEvents_[2][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on MET
	double relMet = met/sumEt;
	hCuts_[3][i]->Fill(relMet);
	if( relMet > cutValue_[3] ) {
	  nCutEvents_[3][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet pt
	hCuts_[4][i]->Fill(jetPt[1]);
	if( jetPt[1] < cutValue_[4] ) {
	  nCutEvents_[4][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet eta
	hCuts_[5][i]->Fill(jetEta[0]);
	hCuts_[5][i]->Fill(jetEta[1]);
	if( std::abs(jetEta[0]) > cutValue_[5] || std::abs(jetEta[1]) > cutValue_[5] ) {
	  nCutEvents_[5][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet Delta phi
	double deltaPhi = TVector2::Phi_mpi_pi(jetPhi[0] - jetPhi[1]);
	hCuts_[6][i]->Fill(deltaPhi);
	if( std::abs(deltaPhi) < cutValue_[6] ) {
	  nCutEvents_[6][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet EMF
	hCuts_[7][i]->Fill(jetEMF[0]);
	hCuts_[7][i]->Fill(jetEMF[1]);
	if( jetEMF[0] < cutValue_[7] || jetEMF[1] < cutValue_[7] ) {
	  nCutEvents_[7][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet n90Hits
	hCuts_[8][i]->Fill(jetN90Hits[0]);
	hCuts_[8][i]->Fill(jetN90Hits[1]);
	if( jetN90Hits[0] < cutValue_[8] || jetN90Hits[1] < cutValue_[8] ) {
	  nCutEvents_[8][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet FRBX
	hCuts_[9][i]->Fill(jetFRBX[0]);
	hCuts_[9][i]->Fill(jetFRBX[1]);
	if( jetFRBX[0] > cutValue_[9] || jetFRBX[1] > cutValue_[9] ) {
	  nCutEvents_[9][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;

	// Cut on jet FHPD
	hCuts_[10][i]->Fill(jetFHPD[0]);
	hCuts_[10][i]->Fill(jetFHPD[1]);
	if( jetFHPD[0] > cutValue_[10] || jetFHPD[1] > cutValue_[10] ) {
	  nCutEvents_[10][i]++;
	  isGood = false;
	}
	if( !isGood ) continue;


	// Fill distributions of selected events
	nGoodEvents_[i]++;
	double diff = jetPt[0]-jetPt[1];
	if( rand_->Uniform() > 0.5 ) diff = jetPt[1]-jetPt[0];
	hPtAsym_[i]->Fill(diff/(jetPt[0]+jetPt[1]));
	if( deltaPhi < 0 ) deltaPhi += 2*M_PI;
	hDeltaPhi_[i]->Fill(deltaPhi);
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

    rand_ = new TRandom3(0);
    treeName_ = treeName;
    nSamples_ = nSamples;

    chain_ = std::vector<TChain*>(nSamples_);
    drawOption_ = std::vector<TString>(nSamples_);
    markerStyle_ = std::vector<int>(nSamples_);
    color_ = std::vector<int>(nSamples_);

    for(int i = 0; i < nSamples_; i++) {
      chain_[i] = new TChain(treeName_);
      
      drawOption_[i] = "HIST";
      markerStyle_[i] = 20;
      color_[i] = 1;
    }

    // Cuts
    nCuts_ = 11;
    cutValue_ = std::vector<double>(nCuts_);
    cutValue_[0] = 2;     //MinNJets
    cutValue_[1] = 1;     //MinVtxNTracks
    cutValue_[2] = 20.;   //MaxVtxPosZ
    cutValue_[3] = 0.5;   //MaxRelMet
    cutValue_[4] = 4.;   //MinJetPt
    cutValue_[5] = 3.;    //MaxJetEta
    cutValue_[6] = 2.1;   //MinDeltaPhi
    cutValue_[7] = 0.01;  //MinEMF
    cutValue_[8] = 2;     //MinN90Hits
    cutValue_[9] = 0.98;  //MaxFHPD
    cutValue_[10] = 0.98; //MaxFRBX

    nEvents_ = std::vector<int>(nSamples_);
    nGoodEvents_ = std::vector<int>(nSamples_);
    nCutLumiBlock_ = std::vector<int>(nSamples_);
    for(int i = 0; i < nSamples_; i++) {
      nEvents_[i] = 0;
      nGoodEvents_[i] = 0;
      nCutLumiBlock_[i] = 0;
    }

    nCutEvents_ = std::vector< std::vector<int> >(nCuts_);
    for(int c = 0; c < nCuts_; c++) {
      nCutEvents_.at(c) = std::vector<int>(nSamples_,0);
    }

    // N-1 distributions
    hCuts_ = std::vector< std::vector<TH1D*> >(nCuts_);
    for(int c = 0; c < nCuts_; c++) {
      hCuts_.at(c) = std::vector<TH1D*>(nSamples_);
      for(int i = 0; i < nSamples_; i++) {
	TString name = "hCuts";
	name += c;
	name += "_";
	name += i;
	TString title = "Cut (";
	title += c;
	title += "): n - 1 plot";
	if( c == 0 ) {
	  hCuts_[c][i] = new TH1D(name,title,30,0,30);
	  hCuts_[c][i]->GetXaxis()->SetTitle("Number of jets");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of events");
	} else if( c == 1 ) {
	  hCuts_[c][i] = new TH1D(name,title,50,0,50);
	  hCuts_[c][i]->GetXaxis()->SetTitle("Number of primary vertex tracks");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of events");
	} else if( c == 2 ) {
	  hCuts_[c][i] = new TH1D(name,title,30,-30,30);
	  hCuts_[c][i]->GetXaxis()->SetTitle("Primary vertex z-position (cm)");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of events");
	} else if( c == 3 ) {
	  hCuts_[c][i] = new TH1D(name,title,50,0,1);
	  hCuts_[c][i]->GetXaxis()->SetTitle("MET / sumET");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of events");
	} else if( c == 4 ) {
	  hCuts_[c][i] = new TH1D(name,title,30,0,30);
	  hCuts_[c][i]->GetXaxis()->SetTitle("p^{jet2}_{T} (GeV)");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of jets");
	} else if( c == 5 ) {
	  hCuts_[c][i] = new TH1D(name,title,20,-5.2,5.2);
	  hCuts_[c][i]->GetXaxis()->SetTitle("#eta");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of jets");
	} else if( c == 6 ) {
	  hCuts_[c][i] = new TH1D(name,title,20,-3.2,3.2);
	  hCuts_[c][i]->GetXaxis()->SetTitle("#Delta #phi");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of jets");
	} else if( c == 7 ) {
	  hCuts_[c][i] = new TH1D(name,title,20,0,1);
	  hCuts_[c][i]->GetXaxis()->SetTitle("f_{EM}");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of jets");
	} else if( c == 8 ) {
	  hCuts_[c][i] = new TH1D(name,title,25,0,50);
	  hCuts_[c][i]->GetXaxis()->SetTitle("n90Hits");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of jets");
	} else if( c == 9 ) {
	  hCuts_[c][i] = new TH1D(name,title,20,0,1);
	  hCuts_[c][i]->GetXaxis()->SetTitle("f_{HPD}");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of jets");
	} else if( c == 10 ) {
	  hCuts_[c][i] = new TH1D(name,title,20,0,1);
	  hCuts_[c][i]->GetXaxis()->SetTitle("f_{RBX}");
	  hCuts_[c][i]->GetYaxis()->SetTitle("Number of jets");
	}

	hCuts_[c][i]->SetMarkerStyle(markerStyle_[i]);
      }
    }

    // Histos
    hPt_ = std::vector<TH1D*>(nSamples_);
    hPtAsym_ = std::vector<TH1D*>(nSamples_);
    hPtCorr_ = std::vector<TH1D*>(nSamples_);
    hEta_ = std::vector<TH1D*>(nSamples_);
    hPhi_ = std::vector<TH1D*>(nSamples_);
    hDeltaPhi_ = std::vector<TH1D*>(nSamples_);
    for(int i = 0; i < nSamples_; i++) {
      chain_[i] = new TChain(treeName_);
      
      TString name;
      name = "hPt";
      name += i;
      hPt_[i] = new TH1D(name,"Selected dijet events",15,0,60);
      hPt_[i]->GetXaxis()->SetTitle("p_{T} (GeV)");
      hPt_[i]->GetYaxis()->SetTitle("Number of jets");
      hPt_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hPtAsym";
      name += i;
      hPtAsym_[i] = new TH1D(name,"Selected dijet events",15,-1.5,1.5);
      hPtAsym_[i]->GetXaxis()->SetTitle("p_{T} asymmetry");
      hPtAsym_[i]->GetYaxis()->SetTitle("Number of dijet events");
      hPtAsym_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hPtCorr";
      name += i;
      hPtCorr_[i] = new TH1D(name,"Selected dijet events",15,0,60);
      hPtCorr_[i]->GetXaxis()->SetTitle("L2L3 corrected p_{T} (GeV)");
      hPtCorr_[i]->GetYaxis()->SetTitle("Number of jets");
      hPtCorr_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hEta";
      name += i;
      hEta_[i] = new TH1D(name,"Selected dijet events",15,-5.,5.);
      hEta_[i]->GetXaxis()->SetTitle("#eta");
      hEta_[i]->GetYaxis()->SetTitle("Number of jets");
      hEta_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hPhi";
      name += i;
      hPhi_[i] = new TH1D(name,"Selected dijet events",15,-3.2,3.2);
      hPhi_[i]->GetXaxis()->SetTitle("#phi");
      hPhi_[i]->GetYaxis()->SetTitle("Number of jets");
      hPhi_[i]->SetMarkerStyle(markerStyle_[i]);

      name = "hDeltaPhi";
      name += i;
      hDeltaPhi_[i] = new TH1D(name,"Selected dijet events",15,2,4);
      hDeltaPhi_[i]->GetXaxis()->SetTitle("#Delta#phi");
      hDeltaPhi_[i]->GetYaxis()->SetTitle("Number of events");
      hDeltaPhi_[i]->SetMarkerStyle(markerStyle_[i]);
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


void normaliseDistributions(int normSample) {
  std::cout << "Normalising distributions to sample '" << normSample << "'... " << std::flush;

  for(int i = 0; i < nSamples_; i++) {
    if( i == normSample ) continue;

    for(int k = 0; k < 6; k++) {
      TH1D *h0 = hPt_.at(normSample);
      TH1D *h1 = hPt_.at(i);
      if( k == 1 ) {
	h0 = hPtAsym_.at(normSample);
	h1 = hPtAsym_.at(i);
      } else if( k == 2 ) {
	h0 = hPtCorr_.at(normSample);
	h1 = hPtCorr_.at(i);
      } else if( k == 3 ) {
	h0 = hEta_.at(normSample);
	h1 = hEta_.at(i);
      } else if( k == 4 ) {
	h0 = hPhi_.at(normSample);
	h1 = hPhi_.at(i);
      } else if( k == 5 ) {
	h0 = hDeltaPhi_.at(normSample);
	h1 = hDeltaPhi_.at(i);
      }

      double norm = h0->Integral();
      if( h1->Integral() ) {
	norm /= h1->Integral();
	h1->Scale(norm);
      }
    }

    for(int c = 0; c < nCuts_; c++) {
      TH1D *h0 = hCuts_.at(c).at(normSample);
      TH1D *h1 = hCuts_.at(c).at(i);
      
      double norm = h0->Integral();
      if( h1->Integral() ) {
	norm /= h1->Integral();
	h1->Scale(norm);
      }
    }
  }
  std::cout << "ok\n";
}


void printCutFlow(int normSample) {
  std::vector<TString> cutLabel(nCuts_);
  cutLabel[0] = "N jets                 >=  ";
  cutLabel[1] = "N vtx tracks           >=  ";
  cutLabel[2] = "Vtx |z|                <   ";
  cutLabel[3] = "MET/SumEt              <   ";
  cutLabel[4] = "Pt jet (2)             >   ";
  cutLabel[5] = "|eta| jet (1,2)        <   ";
  cutLabel[6] = "|DeltaPhi| jet (1,2)   >   ";
  cutLabel[7] = "fEMF jet (1,2)         >   ";
  cutLabel[8] = "n90Hits jet (1,2)      >=  ";
  cutLabel[9] = "fHPD jet (1,2)         <   ";
  cutLabel[10] = "fRBX jet (1,2)        <   ";

  std::cout << "\n\n                                 \tData\tMC\n";

  std::vector<int> nEvts = nEvents_;
  std::cout << "Total                            \t";
  for(int n = 0; n < nSamples_; n++) {
    if( n != normSample ) std::cout << "(";
    std::cout << nEvts[n];
    if( n != normSample ) std::cout << ")";
    std::cout << "\t";
  }
  std::cout << std::endl;
  std::cout << "Good luminosity block               \t";
  for(int n = 0; n < nSamples_; n++) {
    if( n != normSample ) std::cout << "(";
    std::cout << (nEvts[n] -= nCutLumiBlock_[n]);
    if( n != normSample ) std::cout << ")";
    std::cout << "\t";
  }
  std::cout << std::endl;
  for(int c = 0; c < nCuts_; c++) {
    std::cout << c << ": " << cutLabel[c] << cutValue_[c] << "   \t" << std::flush;
    for(int n = 0; n < nSamples_; n++) {
    if( n != normSample ) std::cout << "(";
    std::cout << (nEvts[n] -= nCutEvents_[c][n]);
    if( n != normSample ) std::cout << ")";
    std::cout << "\t";
    }
    std::cout << std::endl;
  }
}
