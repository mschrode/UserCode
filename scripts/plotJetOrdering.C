#include <fstream>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2D.h"
#include "TString.h"
#include "TStyle.h"
#include "TVector2.h"

#include "../util/utils.h"
#include "../util/StyleSettings.h"



// --------------------------------------------------
class JetIndex {
public:
  JetIndex(unsigned int idx, double pt) : idx_(idx), pt_(pt) {};
  const unsigned int idx_;
  const double pt_;
  // For sorting jets in pt
  static bool ptGreaterThan(const JetIndex *idx1, const JetIndex *idx2) {
    // check for 0
    if(idx1 == 0) {
      return idx2 != 0;
    } else if (idx2 == 0) {
      return false;
    } else {
      return idx1->pt_ > idx2->pt_;
    }
  }
};



// --------------------------------------------------
TChain *createTChain(const TString &fileListName) {
  TChain* chain = new TChain("DiJetTree"); 
  std::ifstream filelist;
  filelist.open(fileListName);
  int nOpenedFiles = 0;
  if( filelist.is_open() ) {
    TString name = "";
    while( !filelist.eof() ) {
      filelist >> name;
      if( filelist.eof() ) break;
      chain->Add(name);
      nOpenedFiles++;
    }
  } else {
    std::cerr << "ERROR opening file '" << fileListName << "'\n";
    exit(1);
  }
  filelist.close();
  return chain;
}



// --------------------------------------------------
void plotJetOrdering(int nMaxEvts = -1) {

  util::StyleSettings::presentation();
  gStyle->SetPadTopMargin(0.17);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetTitleYOffset(1.);


  const TString fileName = "input/Kalibri_Calo_Fall10";
  const double minEta = 0.;
  const double maxEta = 1.1;
  const double minDeltaPhi = 2.7;
  const unsigned int jIdxSize = 10;
  const int maxNJet = 50;

  std::vector<double> ptBinEdges;
  ptBinEdges.push_back(90.);
  ptBinEdges.push_back(120.);
  ptBinEdges.push_back(200.);
  ptBinEdges.push_back(500.);
  ptBinEdges.push_back(1000.);


  // Set up histograms
  std::cout << "Setting up histograms" << std::endl;

  std::vector<TH2*> h2Idx(ptBinEdges.size()-1);
  for(size_t i = 0; i < h2Idx.size(); ++i) {
    h2Idx[i] = new TH2D("h2Idx"+util::toTString(i),util::toTString(ptBinEdges[i])+" < p^{ave}_{T} < "+util::toTString(ptBinEdges[i+1])+" GeV;Idx(Uncorr);Idx(Corr)",jIdxSize,-0.5,jIdxSize-0.5,jIdxSize,-0.5,jIdxSize-0.5);
    h2Idx[i]->GetXaxis()->SetNdivisions(10);
    h2Idx[i]->GetYaxis()->SetNdivisions(10);
  }


  // Read variables
  std::cout << "Setting up chain" << std::endl;

  std::vector<JetIndex*> jIdx;
  jIdx.resize(jIdxSize,0);

  float weight = 1.;
  int nObjJet = 0;
  float jetPt[maxNJet];
  float jetEta[maxNJet];
  float jetPhi[maxNJet];
  bool jetID[maxNJet];
  float jetCorrL2L3[maxNJet];

  TChain *chain = createTChain(fileName);
  chain->SetBranchAddress("Weight",&weight);
  chain->SetBranchAddress("NobjJet",&nObjJet);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetEta",jetEta);
  chain->SetBranchAddress("JetPhi",jetPhi);
  chain->SetBranchAddress("JetIDLoose",jetID);
  chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);

  // Loop over tree entries and fill histograms
  std::cout << "Reading events from chain" << std::endl;
  int nEntries = chain->GetEntries();
  if( nMaxEvts > 0 && nEntries > nMaxEvts ) nEntries = nMaxEvts;
  for(int n = 0; n < nEntries; ++n) {
    if( n%10000 == 0 ) std::cout << " Entry " << n << std::endl;
    chain->GetEntry(n);

    if( nObjJet > maxNJet ) {
      std::cerr << "WARNING: nObjJet = " << nObjJet << " > " << maxNJet << ". Skipping event.\n";
      continue;
    }

    if( nObjJet < 3 ) continue;
    

    // Sort by corrected pt
    unsigned int nJets = nObjJet;
    if( jIdxSize < nJets ) nJets = jIdxSize;
    for(size_t i = 0; i < nJets; ++i) {
      jIdx[i] = new JetIndex(i,jetPt[i]*jetCorrL2L3[i]);
    }
    std::sort(jIdx.begin(),jIdx.begin()+nJets,JetIndex::ptGreaterThan);

    // Compute auxiliary quantities
    double ptAveCorr = 0.5*(jIdx[0]->pt_+jIdx[1]->pt_);

    // Selection
    if( ptAveCorr < ptBinEdges.front() || ptAveCorr > ptBinEdges.back() ) continue;
    if( std::abs(jetEta[jIdx[0]->idx_]) < minEta || std::abs(jetEta[jIdx[1]->idx_]) < minEta ||
	std::abs(jetEta[jIdx[0]->idx_]) > maxEta || std::abs(jetEta[jIdx[1]->idx_]) > maxEta ) continue;
    if( std::abs(TVector2::Phi_mpi_pi(jetPhi[jIdx[0]->idx_]-jetPhi[jIdx[1]->idx_])) < minDeltaPhi ) continue;
    if( !(jetID[jIdx[0]->idx_] && jetID[jIdx[1]->idx_]) ) continue;

    // Fill histograms
    unsigned int ptAveCorrBin = 0;
    if( util::findBin(ptAveCorr,ptBinEdges,ptAveCorrBin) ) {
      TH2 *h = h2Idx.at(ptAveCorrBin);
      for(unsigned int i = 0; i < nJets; ++i) {
	h->Fill(i,jIdx[i]->idx_);
      }
    } else {
      std::cerr << "WARNING: bin out-of-range" << std::endl;
    }
  } // End of loop over entries


  // Plot histograms
  for(size_t i = 0; i < h2Idx.size(); ++i) {
    TCanvas *can = new TCanvas("can"+util::toTString(i),"Bin "+util::toTString(i),500,500);
    can->cd();
    h2Idx[i]->Draw("COLZ");
    can->SetLogz();
    can->SaveAs("JetIdxCorrVsUncorr_PtBin"+util::toTString(i)+".eps","eps");
  }
}

