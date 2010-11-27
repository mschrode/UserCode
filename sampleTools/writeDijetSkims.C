// $Id: writeDijetSkims.C,v 1.1 2010/11/26 22:42:31 mschrode Exp $

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"
#include "TVector2.h"

#include "BinningAdmin.h"



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
  std::cout << "Getting trees from ntuples in '" << fileListName << "'" << std::endl;
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
  std::cout << "  Done. Added " << nOpenedFiles << " trees." << std::endl;

  return chain;
}



// --------------------------------------------------
void writeDijetSkims() {


  // ++++ Set parameters +++++++++++++++++++++++++++++++++++++++

  const int nEvts = -10000;
  const TString config = "input/BinningAdmin.cfg";
  const TString inFileListName = "input/Kalibri_PF_Fall10";
  //  const TString inFileListName = "input/Kalibri_Calo_0200";
  const bool isData = false;
  const unsigned int maxHltThres = 140;
  const double minDeltaPhi = 2.7;



  // ++++ Checks and follow-up parameter +++++++++++++++++++++++

  // Prepare name of output files  
  TString outFilePrefix = "~/lustre/KalibriDiJetSkims/KalibriDiJetSkim";

  if( inFileListName.Contains("Calo") ) outFilePrefix += "_Calo";
  else if( inFileListName.Contains("PF") ) outFilePrefix += "_PF";  

  if( isData ) outFilePrefix += "_Data";
  else if( inFileListName.Contains("Fall10") ) outFilePrefix += "_MCFall10";  
  else outFilePrefix += "_MC";  


  sampleTools::BinningAdmin binAdmin(config);

  const TString hlt = isData ? binAdmin.triggerName(maxHltThres) : "none";

  binAdmin.printBinning();
  if( isData ) binAdmin.print(hlt);
  
  unsigned int nMaxNJet = 0;
  unsigned int nDijets = 0;
  unsigned int nHlt = 0;
  unsigned int nDeltaPhi = 0;
  unsigned int nJetID = 0;
  unsigned int nEta = 0;
  unsigned int nPtAve = 0;

  


  // ++++ Prepare trees and files ++++++++++++++++++++++++++++++

  // Open Kalibri ntuples
  TChain *oldChain = createTChain(inFileListName);

  //TChain *oldChain = new TChain("DiJetTree");
  //oldChain->Add("QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6job_0_ak5Calo.root");
  
  // Deactivate branches not needed
  oldChain->SetBranchStatus("Track*",0);
  oldChain->SetBranchStatus("NobjTrack",0);
  oldChain->SetBranchStatus("Tow*",0);
  oldChain->SetBranchStatus("GenPart*",0);
  oldChain->SetBranchStatus("GenPartId*",0);
  oldChain->SetBranchStatus("Vtx*",0);
  oldChain->SetBranchStatus("VtxN",1);
  oldChain->SetBranchStatus("Met*",0);
  oldChain->SetBranchStatus("Mu*",0);

  // Get elements needed for preselection
  const int maxNJet = 50;
  int nObjJet = 0;
  float jetPt[maxNJet];
  float jetEta[maxNJet];
  float jetPhi[maxNJet];
  float jetCorrL2L3[maxNJet];
  bool jetID[maxNJet];
  bool hltDiJetAve50U = false;
  bool hltDiJetAve70U = false;
  bool hltDiJetAve100U = false;
  bool hltDiJetAve140U = false;

  // Get other elements written to tree
  Float_t         JetEt[maxNJet];   //[NobjJet]
  Float_t         JetE[maxNJet];   //[NobjJet]
  Int_t           JetN90Hits[maxNJet];   //[NobjJet]
  Float_t         JetHad[maxNJet];   //[NobjJet]
  Float_t         JetEMF[maxNJet];   //[NobjJet]
  Float_t         JetFHPD[maxNJet];   //[NobjJet]
  Float_t         JetFRBX[maxNJet];   //[NobjJet]
  Bool_t          JetIDTight[maxNJet];   //[NobjJet]
  Float_t         JetEtWeightedSigmaPhi[maxNJet];   //[NobjJet]
  Float_t         JetEtWeightedSigmaEta[maxNJet];   //[NobjJet]
  Float_t         JetCorrZSP[maxNJet];   //[NobjJet]
  Float_t         JetCorrL2[maxNJet];   //[NobjJet]
  Float_t         JetCorrL3[maxNJet];   //[NobjJet]
  Float_t         JetCorrJPT[maxNJet];   //[NobjJet]
  Float_t         JetCorrL2L3JPT[maxNJet];   //[NobjJet]
  Float_t         JetCorrL4JW[maxNJet];   //[NobjJet]
  Int_t           JetIEta[maxNJet];   //[NobjJet]
  Int_t           JetIPhi[maxNJet];   //[NobjJet]
  Float_t         JetGenJetDeltaR[maxNJet];   //[NobjJet]
  Float_t         GenJetPt[maxNJet];   //[NobjJet]
  Float_t         GenJetPhi[maxNJet];   //[NobjJet]
  Float_t         GenJetEta[maxNJet];   //[NobjJet]
  Float_t         GenJetEt[maxNJet];   //[NobjJet]
  Float_t         GenJetE[maxNJet];   //[NobjJet]
  Int_t           NobjGenJet;
  Float_t         GenJetColPt[maxNJet];   //[NobjGenJet]
  Float_t         GenJetColPhi[maxNJet];   //[NobjGenJet]
  Float_t         GenJetColEta[maxNJet];   //[NobjGenJet]
  Float_t         GenJetColEt[maxNJet];   //[NobjGenJet]
  Float_t         GenJetColE[maxNJet];   //[NobjGenJet]
  Float_t         GenJetColEmE[maxNJet];   //[NobjGenJet]
  Float_t         GenJetColHadE[maxNJet];   //[NobjGenJet]
  Float_t         GenJetColInvE[maxNJet];   //[NobjGenJet]
  Float_t         GenJetColAuxE[maxNJet];   //[NobjGenJet]
  Int_t           GenJetColJetIdx[maxNJet];   //[NobjGenJet]


  // Set branch addresses
  oldChain->SetBranchAddress("NobjJet",&nObjJet);
  oldChain->SetBranchAddress("JetPt",jetPt);
  oldChain->SetBranchAddress("JetEta",jetEta);
  oldChain->SetBranchAddress("JetPhi",jetPhi);
  oldChain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  oldChain->SetBranchAddress("JetIDLoose",jetID);
  oldChain->SetBranchAddress("HltDiJetAve50U",&hltDiJetAve50U);
  oldChain->SetBranchAddress("HltDiJetAve70U",&hltDiJetAve70U);
  oldChain->SetBranchAddress("HltDiJetAve100U",&hltDiJetAve100U);
  oldChain->SetBranchAddress("HltDiJetAve140U",&hltDiJetAve140U);

  oldChain->SetBranchAddress("JetEt", JetEt);
  oldChain->SetBranchAddress("JetE", JetE);
  oldChain->SetBranchAddress("JetN90Hits", JetN90Hits);
  oldChain->SetBranchAddress("JetHad", JetHad);
  oldChain->SetBranchAddress("JetEMF", JetEMF);
  oldChain->SetBranchAddress("JetFHPD", JetFHPD);
  oldChain->SetBranchAddress("JetFRBX", JetFRBX);
  oldChain->SetBranchAddress("JetIDTight", JetIDTight);
  oldChain->SetBranchAddress("JetEtWeightedSigmaPhi", JetEtWeightedSigmaPhi);
  oldChain->SetBranchAddress("JetEtWeightedSigmaEta", JetEtWeightedSigmaEta);
  oldChain->SetBranchAddress("JetCorrZSP", JetCorrZSP);
  oldChain->SetBranchAddress("JetCorrL2", JetCorrL2);
  oldChain->SetBranchAddress("JetCorrL3", JetCorrL3);
  oldChain->SetBranchAddress("JetCorrJPT", JetCorrJPT);
  oldChain->SetBranchAddress("JetCorrL2L3JPT", JetCorrL2L3JPT);
  oldChain->SetBranchAddress("JetCorrL4JW", JetCorrL4JW);
  oldChain->SetBranchAddress("JetIEta", JetIEta);
  oldChain->SetBranchAddress("JetIPhi", JetIPhi);
  oldChain->SetBranchAddress("JetGenJetDeltaR", JetGenJetDeltaR);
  oldChain->SetBranchAddress("GenJetPt", GenJetPt);
  oldChain->SetBranchAddress("GenJetPhi", GenJetPhi);
  oldChain->SetBranchAddress("GenJetEta", GenJetEta);
  oldChain->SetBranchAddress("GenJetEt", GenJetEt);
  oldChain->SetBranchAddress("GenJetE", GenJetE);
  oldChain->SetBranchAddress("NobjGenJet", &NobjGenJet);
  oldChain->SetBranchAddress("GenJetColPt", GenJetColPt);
  oldChain->SetBranchAddress("GenJetColPhi", GenJetColPhi);
  oldChain->SetBranchAddress("GenJetColEta", GenJetColEta);
  oldChain->SetBranchAddress("GenJetColEt", GenJetColEt);
  oldChain->SetBranchAddress("GenJetColE", GenJetColE);
  oldChain->SetBranchAddress("GenJetColEmE", GenJetColEmE);
  oldChain->SetBranchAddress("GenJetColHadE", GenJetColHadE);
  oldChain->SetBranchAddress("GenJetColInvE", GenJetColInvE);
  oldChain->SetBranchAddress("GenJetColAuxE", GenJetColAuxE);
  oldChain->SetBranchAddress("GenJetColJetIdx", GenJetColJetIdx);
  

  //Create a new file + a clone of old tree in new file per eta and pt bin
  std::vector< std::vector<TFile*> > newFiles(binAdmin.nEtaBins());
  std::vector< std::vector<TTree*> > newTrees(binAdmin.nEtaBins());
  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    newFiles[etaBin] = std::vector<TFile*>(binAdmin.nPtBins(hlt,etaBin));
    newTrees[etaBin] = std::vector<TTree*>(binAdmin.nPtBins(hlt,etaBin));
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
      TString name = outFilePrefix+"_Eta";
      name += etaBin;
      name += "_Pt";
      name += (binAdmin.hltMinPtBin(hlt,etaBin)+ptBin);
      name += ".root";
      newFiles[etaBin][ptBin] = new TFile(name,"RECREATE");
      newTrees[etaBin][ptBin] = oldChain->CloneTree(0);
    }
  }

  // Add branch with indices of jets ordered by
  // L2L3 corrected pt to new tree
  int corrJetIdx[maxNJet];
  for(int j = 0; j < maxNJet; ++j) {
    corrJetIdx[j] = -1;
  }
  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
      newTrees[etaBin][ptBin]->Branch("L2L3CorrJetColJetIdx",corrJetIdx,"L2L3CorrJetColJetIdx[NobjJet]/I");
    }
  }

  // Struct for jet ordering
  std::vector<JetIndex*> jIdx(maxNJet);




  // ++++ Loop over old tree and select dijets +++++++++++++++++

  int nEntries = oldChain->GetEntries();
  if( nEvts > 0 && nEvts <= nEntries ) nEntries = nEvts;

  for(int i = 0; i < nEntries; ++i) {
    oldChain->GetEntry(i);

    if( nObjJet > maxNJet ) {
      std::cerr << "WARNING: nObjJet = " << nObjJet << " > " << maxNJet << ". Skipping event.\n";
      ++nMaxNJet;
      continue;
    }

    if( nObjJet < 2 ) {
      ++nDijets;
      continue;
    }

    if( isData ) {
      // HLT cuts
      if( maxHltThres == 50 && !hltDiJetAve50U ) 
 	continue;
      else if( maxHltThres == 70 && !(hltDiJetAve50U || hltDiJetAve70U) ) 
 	continue;
      else if( maxHltThres == 100 && !(hltDiJetAve50U || hltDiJetAve70U || hltDiJetAve100U) ) 
 	continue;
      else if( maxHltThres == 140 && !(hltDiJetAve50U || hltDiJetAve70U || hltDiJetAve100U || hltDiJetAve140U) ) 
 	continue;
    }
  
    // Order L2L3 corrected jets
    for(int j = 0; j < nObjJet; ++j) {
      jIdx[j] = new JetIndex(j,jetPt[j]*jetCorrL2L3[j]);
    }
    std::sort(jIdx.begin(),jIdx.begin()+nObjJet,JetIndex::ptGreaterThan);
    for(int j = 0; j < nObjJet; ++j) {
      corrJetIdx[j] = jIdx[j]->idx_;
    }
    for(int j = 0; j < nObjJet; ++j) {
      delete jIdx[j];
    }
  
    if( std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJetIdx[0]]-jetPhi[corrJetIdx[1]])) < minDeltaPhi ) {
      ++nDeltaPhi;
      continue;
    }
  
    if( !jetID[corrJetIdx[0]] || !jetID[corrJetIdx[1]] ) {
      ++nJetID;
      continue;
    }

    // Find eta and pt bin
    unsigned int etaBin = 1000;
    if( binAdmin.findSameEtaBin(jetEta[corrJetIdx[0]],jetEta[corrJetIdx[1]],etaBin) ) {
      double ptAve = 0.5*(jetCorrL2L3[corrJetIdx[0]]*jetPt[corrJetIdx[0]]+jetCorrL2L3[corrJetIdx[1]]*jetPt[corrJetIdx[1]]);
      unsigned int ptAveBin = 1000;
      if( binAdmin.findPtBin(hlt,ptAve,etaBin,ptAveBin) ) {
 	ptAveBin -= binAdmin.hltMinPtBin(hlt,etaBin);
 	newTrees[etaBin][ptAveBin]->Fill();
      } else {
 	++nPtAve;
      }
    } else {
      ++nEta;
    }

    if( i%50000 == 0 ) {
      std::cout << "Processed " << i << " events" << std::endl;
      for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
 	for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
 	  newTrees[etaBin][ptBin]->AutoSave();
 	}
      }
    }
  } // End of loop over entries

  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
      newTrees[etaBin][ptBin]->AutoSave();
      //newTrees[etaBin][ptBin]->Print();
    }
  }




  // ++++ Print status ++++++++++++++++++++++++++++++++++++++++++
  std::cout << "Done processing " << nEntries << " events from file list '" << inFileListName << "'" << std::endl;

  std::cout << "Selected " << std::endl;
  std::cout << "  " << (nEntries -= nMaxNJet ) << " events with <= " << maxNJet << " jets " << std::endl;
  std::cout << "  " << (nEntries -= nDijets ) << " events with >= 2 jets " << std::endl;
  if( isData ) std::cout << "  " << (nEntries -= nHlt ) << " events passing HLT trigger (max threshold " << maxHltThres << " GeV)" << std::endl;
  std::cout << "  " << (nEntries -= nDeltaPhi ) << " events with |DeltaPhi(1,2)| > " << minDeltaPhi << std::endl;
  std::cout << "  " << (nEntries -= nJetID ) << " events with Jet(1,2) passing loose JetID cuts" << std::endl;
  std::cout << "  " << (nEntries -= nPtAve ) << " events with PtAve within binning" << std::endl;
  std::cout << "  " << (nEntries -= nEta ) << " events with Eta(1,2) within binning" << std::endl;

  std::cout << "Wrote " << nEntries << " events to files" << std::endl;
  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
      std::cout << "  " << newTrees[etaBin][ptBin]->GetEntries() << " events with " << binAdmin.etaMin(etaBin) << " < |eta(1,2)| < " << binAdmin.etaMax(etaBin) << " to file '" << newFiles[etaBin][ptBin]->GetName() << "'" << std::endl;
    }
  }
    


  // ++++ Clean up +++++++++++++++++++++++++++++++++++++++++++++

  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
      newFiles[etaBin][ptBin]->Close();
      delete newFiles[etaBin][ptBin];
    }
  }
  delete oldChain;
}
