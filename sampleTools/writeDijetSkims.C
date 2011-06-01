// $Id: writeDijetSkims.C,v 1.4 2011/05/20 10:00:07 mschrode Exp $
//
// Skim Kalibri ntuples as input for resolution fit.
// At this pre-selection
//  1) data events are selected from fully efficient trigger
//     paths;
//  2) jets are ordered in L1*L2*L3 corrected pt. A new branch
//     'L2L3CorrJetColJetIdx' is included into the ntuple, storing
//     for each jet (ordered in raw pt) the index of this jet
//     in the corresponding collection of corrected jet pt;
//  3) dijet events are selected by requiring |DeltaPhi(1,2)| > 2.7;
//  4) the leading two jets are required to pass the loose jet id.
// Then, events are sorted into bins of eta and ptAve as specified
// by an input config file. Per eta and ptAve bin, a separate file
// is created.


#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TBranch.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TString.h"
#include "TVector2.h"

#include "BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"



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
void writeDijetSkims(bool isData, unsigned int maxHltThres = 0) {


  // ++++ Set parameters +++++++++++++++++++++++++++++++++++++++

  const int nEvts = -100;
  const TString config = "BinningAdmin.cfg";
  const TString inFileListName = "input/Analysis2011/Kalibri_MCSummer11_QCDFlat_PUS3_L1FastJet_AK5Calo";
  //const bool isData = true;
  //  const unsigned int maxHltThres = 175;
  const double minDeltaPhi = 2.7;
  const unsigned int minRunNumber = 163337;



  // ++++ Checks and follow-up parameters ++++++++++++++++++++++

  // Prepare name of output files  
  TString outFilePrefix = "~/lustre/Analysis2011/KalibriDiJetSkims_QCD_Pt-15to3000_TuneD6T_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/KalibriSkim";
  if( inFileListName.Contains("Calo") ) outFilePrefix += "_Calo";
  else if( inFileListName.Contains("PF") ) outFilePrefix += "_PF";  
  else if( inFileListName.Contains("JPT") ) outFilePrefix += "_JPT";  

  if( inFileListName.Contains("L1Offset") ) outFilePrefix += "_L1Offset";
  else if( inFileListName.Contains("L1FastJet") ) outFilePrefix += "_L1FastJet";

  if( isData ) outFilePrefix += "_Data";
  else if( inFileListName.Contains("Fall10") ) outFilePrefix += "_MCFall10";  
  else if( inFileListName.Contains("Spring11") ) outFilePrefix += "_MCSpring11";  
  else if( inFileListName.Contains("Summer11") ) outFilePrefix += "_MCSummer11";  
  else outFilePrefix += "_MC";  


  sampleTools::BinningAdmin binAdmin(config);

  const TString hlt = isData ? binAdmin.triggerName(maxHltThres) : "none";

  binAdmin.printBinning();
  if( isData ) binAdmin.print(hlt);
  
  unsigned int nNewTrig = 0;
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
  float jetCorrL1[maxNJet];
  float jetCorrL2L3[maxNJet];	// Should contain L2*L3*LResidual, starting from L1 corrected pt
  bool jetID[maxNJet];
  bool hlt30 = false;
  bool hlt60 = false;
  bool hlt80 = false;
  bool hlt110 = false;
  bool hlt150 = false;
  bool hlt190 = false;
  bool hlt240 = false;
  bool hlt300 = false;
  bool hlt370 = false;

  // Get other elements written to tree
  UInt_t          RunNumber;
  UInt_t          LumiBlockNumber;
  UInt_t          EventNumber;
  Int_t           NobjTow;
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
  Int_t           VtxN = 0;
  Int_t           PUMCNumVtx = 0;


  // Set branch addresses
  oldChain->SetBranchAddress("RunNumber", &RunNumber);
  oldChain->SetBranchAddress("LumiBlockNumber", &LumiBlockNumber);
  oldChain->SetBranchAddress("EventNumber", &EventNumber);
  oldChain->SetBranchAddress("NobjTow",&NobjTow);
  oldChain->SetBranchAddress("NobjJet",&nObjJet);
  oldChain->SetBranchAddress("JetPt",jetPt);
  oldChain->SetBranchAddress("JetEta",jetEta);
  oldChain->SetBranchAddress("JetPhi",jetPhi);
  oldChain->SetBranchAddress("JetCorrL1",jetCorrL1);
  oldChain->SetBranchAddress("JetCorrL2",JetCorrL2);
  oldChain->SetBranchAddress("JetCorrL3",JetCorrL3);
  oldChain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  oldChain->SetBranchAddress("HltDiJetAve30",&hlt30);
  oldChain->SetBranchAddress("HltDiJetAve60",&hlt60);
  oldChain->SetBranchAddress("HltDiJetAve80",&hlt80);
  oldChain->SetBranchAddress("HltDiJetAve110",&hlt110);
  oldChain->SetBranchAddress("HltDiJetAve150",&hlt150);
  oldChain->SetBranchAddress("HltDiJetAve190",&hlt190);
  oldChain->SetBranchAddress("HltDiJetAve240",&hlt240);
  oldChain->SetBranchAddress("HltDiJetAve300",&hlt300);
  oldChain->SetBranchAddress("HltDiJetAve370",&hlt370);
  oldChain->SetBranchAddress("JetEt", JetEt);
  oldChain->SetBranchAddress("JetE", JetE);
  oldChain->SetBranchAddress("JetN90Hits", JetN90Hits);
  oldChain->SetBranchAddress("JetHad", JetHad);
  oldChain->SetBranchAddress("JetEMF", JetEMF);
  oldChain->SetBranchAddress("JetFHPD", JetFHPD);
  oldChain->SetBranchAddress("JetFRBX", JetFRBX);
  oldChain->SetBranchAddress("JetIDLoose",jetID);
  oldChain->SetBranchAddress("JetIDTight", JetIDTight);
  oldChain->SetBranchAddress("JetEtWeightedSigmaPhi", JetEtWeightedSigmaPhi);
  oldChain->SetBranchAddress("JetEtWeightedSigmaEta", JetEtWeightedSigmaEta);
  oldChain->SetBranchAddress("JetCorrZSP", JetCorrZSP);
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
  oldChain->SetBranchAddress("VtxN",&VtxN);
  oldChain->SetBranchAddress("PUMCNumVtx",&PUMCNumVtx);
  

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
  // corrected pt to new tree
  int corrJetIdx[maxNJet];
  for(int j = 0; j < maxNJet; ++j) {
    corrJetIdx[j] = -1;
  }
  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
      newTrees[etaBin][ptBin]->Branch("L2L3CorrJetColJetIdx",corrJetIdx,"L2L3CorrJetColJetIdx[NobjJet]/I");
//       newTrees[etaBin][ptBin]->Branch("HltDiJetAve30",&hlt30,"HltDiJetAve30/O");
//       newTrees[etaBin][ptBin]->Branch("HltDiJetAve60",&hlt60,"HltDiJetAve60/O");
//       newTrees[etaBin][ptBin]->Branch("HltDiJetAve80",&hlt80,"HltDiJetAve80/O");
//       newTrees[etaBin][ptBin]->Branch("HltDiJetAve110",&hlt110,"HltDiJetAve110/O");
//       newTrees[etaBin][ptBin]->Branch("HltDiJetAve150",&hlt150,"HltDiJetAve150/O");
//       newTrees[etaBin][ptBin]->Branch("HltDiJetAve190",&hlt190,"HltDiJetAve190/O");
//       newTrees[etaBin][ptBin]->Branch("HltDiJetAve240",&hlt240,"HltDiJetAve240/O");
//       newTrees[etaBin][ptBin]->Branch("HltDiJetAve300",&hlt300,"HltDiJetAve300/O");
//       newTrees[etaBin][ptBin]->Branch("HltDiJetAve370",&hlt370,"HltDiJetAve370/O");
    }
  }

  // Container for jet ordering
  util::JetIndexCol corrJets;




  // ++++ Loop over old tree and select dijets +++++++++++++++++

  int nEntries = oldChain->GetEntries();
  if( nEvts > 0 && nEvts <= nEntries ) nEntries = nEvts;

  for(int i = 0; i < nEntries; ++i) {
    if( i%50000 == 0 ) {
      std::cout << "Processed " << i << " events" << std::endl;
    }

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
      // Use only 2011 runs with new triggers
      if( RunNumber < minRunNumber ) {
	++nNewTrig;
	continue;
      }

      // HLT cuts
      if( maxHltThres == 30 && !(hlt30) ) {
	++nHlt;
 	continue;
      } else if( maxHltThres == 60 && !(hlt30 || hlt60) ) {
	++nHlt;
 	continue;
      } else if( maxHltThres == 80 && !(hlt30 || hlt60 || hlt80) ) {
	++nHlt;
 	continue;
      } else if( maxHltThres == 110 && !(hlt30 || hlt60 || hlt80 || hlt110) ) {
	++nHlt;
 	continue;
      } else if( maxHltThres == 150 && !(hlt30 || hlt60 || hlt80 || hlt110 || hlt150) ) {
	++nHlt;
 	continue;
      } else if( maxHltThres == 190 && !(hlt30 || hlt60 || hlt80 || hlt110 || hlt150 || hlt190) ) {
	++nHlt;
 	continue;
      } else if( maxHltThres == 240 && !(hlt30 || hlt60 || hlt80 || hlt110 || hlt150 || hlt190 || hlt240) ) {
	++nHlt;
 	continue;
      } else if( maxHltThres == 300 && !(hlt30 || hlt60 || hlt80 || hlt110 || hlt150 || hlt190 || hlt240 || hlt300) ) {
	++nHlt;
 	continue;
      } else if( maxHltThres == 370 && !(hlt30 || hlt60 || hlt80 || hlt110 || hlt150 || hlt190 || hlt240 || hlt300 || hlt370) ) {
	++nHlt;
 	continue;
      }
    }

    // Order corrected jets
    corrJets.clear();
    for(int j = 0; j < nObjJet; ++j) {
      jetCorrL2L3[j] = JetCorrL2[j]*JetCorrL3[j]; // Remove residual correction due to bug in 42

      corrJets.add(j,jetCorrL1[j]*jetCorrL2L3[j]*jetPt[j]);
    }
    corrJets.sort();

    // Set branch corrected jet indices
    for(int j = 0; j < nObjJet; ++j) {
      corrJetIdx[j] = corrJets(j);
    }

    
    // Dijet selection
    if( std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets(0)]-jetPhi[corrJets(1)])) < minDeltaPhi ) {
      ++nDeltaPhi;
      continue;
    }
    else if( !jetID[corrJets(0)] || !jetID[corrJets(1)] ) {
      ++nJetID;
      continue;
    }

    // Find eta and pt bin
    unsigned int etaBin = 1000;
    if( binAdmin.findSameEtaBin(jetEta[corrJets(0)],jetEta[corrJets(1)],etaBin) ) {
      double ptAve = 0.5*(corrJets.pt(0)+corrJets.pt(1));
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
      for(etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
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
  if( isData ) std::cout << "  " << (nEntries -= nNewTrig ) << " events with run number >= " << minRunNumber << std::endl;
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


void writeDijetSkimsData() {
  std::vector<unsigned int> hltThes;
  hltThes.push_back(30);
  hltThes.push_back(60);
  hltThes.push_back(80);
  hltThes.push_back(110);
  hltThes.push_back(150);
  hltThes.push_back(190);
  hltThes.push_back(240);
  hltThes.push_back(300);
  hltThes.push_back(370);
  for(std::vector<unsigned int>::const_iterator it = hltThes.begin();
      it != hltThes.end(); ++it) {
    writeDijetSkims(true,*it);
  }
}
