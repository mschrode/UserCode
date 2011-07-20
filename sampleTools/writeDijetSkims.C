// $Id: writeDijetSkims.C,v 1.7 2011/07/18 09:12:31 mschrode Exp $
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
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TString.h"
#include "TVector2.h"

#include "BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"


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


// Generate weights for Flat10 PU scenario for given
// data PU distribution
// Code from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
// --------------------------------------------------
std::vector<double> generate_flat10_weights(const TH1* data_npu_estimated) {
  // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
  const double npu_probs[25] = {0.0698146584, 0.0698146584, 0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584 /* <-- 10*/,0.0630151648,0.0526654164,0.0402754482,0.0292988928,0.0194384503,0.0122016783,0.007207042,0.004003637,0.0020278322,0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 /* <-- 24 */};
  std::vector<double> result(25);
  double s = 0.0;
  for(int npu=0; npu<25; ++npu) {
    double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
    result[npu] = npu_estimated / npu_probs[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(int npu=0; npu<25; ++npu) {
    result[npu] /= s;
  }
  return result;
}



// --------------------------------------------------
void writeDijetSkims(bool isData, unsigned int maxHltThres = 0) {
  
  
  // ++++ Set parameters +++++++++++++++++++++++++++++++++++++++
  
  std::cout << "Setting parameters" << std::endl;
  
  const TString inFileListName = "input/Analysis2011/Kalibri_MCSummer11_QCDFlat_PythiaZ2_PUS3_L1FastJet_AK5PF";
  const TString outFilePath = "~/lustre/Analysis2011/tmp";
  const TString config = "BinningAdmin.cfg";
  const bool writeAllTrees = false;
  
  const int nEvts = -10000;
  const TString outFilePrefix = outFilePath+"/KalibriSkim";
  const unsigned int minRunNumber = 163337;

  // Event selection cuts
  const double minDeltaPhi = 2.7;
  const double maxDeltaR = inFileListName.Contains("Calo") ? 0.2 : 0.1;	// Different MC truth matching depending on jet type

  // Cut-flow counters
  unsigned int nNewTrig = 0;
  unsigned int nMaxNJet = 0;
  unsigned int nDijets = 0;
  unsigned int nHlt = 0;
  unsigned int nDeltaPhi = 0;
  unsigned int nJetID = 0;
  unsigned int nEta = 0;
  unsigned int nPtAve = 0;
  
  // Event weights
  double wPU = 1.;		// PU weight
  float wSpec = 1.;		// Spectrum weight (for flat samples, taken from input ntuple)
  
  

  // ++++ Prepare objects ++++++++++++++++++++++++++++++++++++++

  std::cout << "Defining binning" << std::endl;

  const sampleTools::BinningAdmin binAdmin(config);
  const TString hlt = isData ? binAdmin.triggerName(maxHltThres) : "none";
  binAdmin.printBinning();
  if( isData ) binAdmin.print(hlt);

  std::vector<double> weightsPU;
  std::vector<TH2*> hRespVsPtGen; // One entry per eta bin
  std::vector<TH2*> hRespVsPtGenLog; // One entry per eta bin
  std::vector<TH2*> hRespVsPtGenLogNoDeltaR;
  std::vector<TH2*> hRespVsPtGenLogNoDeltaPhi;
  std::vector<TH2*> hRespVsPtGenPUReweighted; // One entry per eta bin
  std::vector<TH2*> hRespVsPtGenLogPUReweighted; // One entry per eta bin
  std::vector<TH2*> hRespVsPtGenLogPULow; // One entry per eta bin
  std::vector<TH2*> hRespVsPtGenLogPUHigh; // One entry per eta bin
  std::vector< std::vector<TH2*> > hRespVsPtGenLogVsPtSoftPUReweighted; // [etaBin][pt3Bin]
  std::vector< std::vector<TH1*> > hPtGen; // [etaBin][pt3Bin]
  
  if( !isData ) {
    std::cout << "Preparing PU reweighting" << std::endl;
    
    TH1* hDataPU = util::FileOps::readTH1("~/Kalibri/input/PUDist_Cert_160404-163869_7TeV_May10ReReco.root","pileup");
    weightsPU = generate_flat10_weights(hDataPU);
    
    
    std::cout << "Preparing histograms" << std::endl;
    
    // MC truth
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      TString title = util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin));
      
      // Pt binning as for resolution measurement
      TH2* h = new TH2D("hRespVsPtGen_Eta"+util::toTString(etaBin),title+";p^{gen}_{T} (GeV);Response",
			binAdmin.nPtBins(etaBin),&(binAdmin.ptBinEdges(etaBin).front()),201,0.,2.);
      h->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", |#Delta#Phi|>"+util::toTString(minDeltaPhi));
      h->SetNdivisions(505);
      h->Sumw2();
      hRespVsPtGen.push_back(h);
      hRespVsPtGenPUReweighted.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenPUReweighted_Eta"+util::toTString(etaBin))));
      hRespVsPtGenPUReweighted.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", |#Delta#Phi|>"+util::toTString(minDeltaPhi)+", PU reweighted");
      
      // Equidistant log bins
      std::vector<double> binEdges(35+1);
      util::HistOps::equidistLogBins(binEdges,binEdges.size()-1,4.,2700.);
      h = new TH2D("hRespVsPtGenLog_Eta"+util::toTString(etaBin),title,
		   binEdges.size()-1,&(binEdges.front()),201,0.,2.);
      h->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", |#Delta#Phi|>"+util::toTString(minDeltaPhi));
      h->SetNdivisions(505);
      h->Sumw2();
      hRespVsPtGenLog.push_back(h);
      hRespVsPtGenLogNoDeltaR.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenLogNoDeltaR_Eta"+util::toTString(etaBin))));
      hRespVsPtGenLogNoDeltaR.back()->SetTitle(title);
      title += ", |#DeltaR|<"+util::toTString(maxDeltaR);
      hRespVsPtGenLogNoDeltaPhi.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenLogNoDeltaPhi_Eta"+util::toTString(etaBin))));
      hRespVsPtGenLogNoDeltaPhi.back()->SetTitle(title);
      title += ", |#Delta#Phi|<"+util::toTString(minDeltaPhi);
      hRespVsPtGenLogPUReweighted.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenLogPUReweighted_Eta"+util::toTString(etaBin))));
      hRespVsPtGenLogPUReweighted.back()->SetTitle(title+", PU reweighted");
      hRespVsPtGenLogPULow.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenLogPULow_Eta"+util::toTString(etaBin))));
      hRespVsPtGenLogPULow.back()->SetTitle(title+", low PU");
      hRespVsPtGenLogPUHigh.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenLogPUHigh_Eta"+util::toTString(etaBin))));
      hRespVsPtGenLogPUHigh.back()->SetTitle(title+", high PU");
      
      std::vector<TH2*> vtmp;
      for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
	title = util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin))+", p^{rel}_{T,3} < "+util::toTString(binAdmin.ptSoftMax(ptSoftBin))+", |#DeltaR|<"+util::toTString(maxDeltaR)+", |#Delta#Phi|>"+util::toTString(minDeltaPhi)+", PU reweighted;p^{gen}_{T} (GeV);Response";
	
	h = new TH2D("hRespVsPtGenLogVsPtSoftPUReweighted_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin),title,binEdges.size()-1,&(binEdges.front()),201,0.,2.);
	h->SetNdivisions(505);
	h->Sumw2();
	vtmp.push_back(h);
      }
      hRespVsPtGenLogVsPtSoftPUReweighted.push_back(vtmp);
    }


    // PtGen spectra
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      std::vector<TH1*> vtmp;
      for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
	TString title = util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin))+", p^{rel}_{T,3} < "+util::toTString(binAdmin.ptSoftMax(ptSoftBin))+";p^{gen}_{T} (GeV);Events";

	TH1* h = new TH1D("hPtGen_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin),
			  title,500,0.,2000.);
	h->GetXaxis()->SetNdivisions(505);
	h->SetLineWidth(2);
	h->Sumw2();
	vtmp.push_back(h);
      }
      hPtGen.push_back(vtmp);
    }
  }


  std::cout << "Preparing trees" << std::endl;

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
  Float_t         JetCorrUncert[maxNJet];
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
  Float_t         Weight = 1.;


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
  oldChain->SetBranchAddress("JetCorrUncert",JetCorrUncert);
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
  oldChain->SetBranchAddress("Weight",&Weight);
  

  // Create a new file + a clone of old tree in new file
  // 1) per eta and ptAve bin;
  // 2) per eta and pt bin;
  // 3) per etaGen and ptGenAve bin.

  // Prepare files and tress
  std::vector< std::vector<TFile*> > newFilesEtaPtAve(binAdmin.nEtaBins());
  std::vector< std::vector<TTree*> > newTreesEtaPtAve(binAdmin.nEtaBins());
  std::vector< std::vector<TFile*> > newFilesEtaPt(binAdmin.nEtaBins());
  std::vector< std::vector<TTree*> > newTreesEtaPt(binAdmin.nEtaBins());
  std::vector< std::vector<TFile*> > newFilesEtaGenPtGenAve(binAdmin.nEtaBins());
  std::vector< std::vector<TTree*> > newTreesEtaGenPtGenAve(binAdmin.nEtaBins());

  // Prepare name of output files  
  TString outFileId = "";
  if( inFileListName.Contains("Calo") ) outFileId += "_Calo";
  else if( inFileListName.Contains("PF") ) outFileId += "_PF";  
  else if( inFileListName.Contains("JPT") ) outFileId += "_JPT";  

  if( inFileListName.Contains("L1Offset") ) outFileId += "_L1Offset";
  else if( inFileListName.Contains("L1FastJet") ) outFileId += "_L1FastJet";

  if( isData ) outFileId += "_Data";
  else if( inFileListName.Contains("Fall10") ) outFileId += "_MCFall10";  
  else if( inFileListName.Contains("Spring11") ) outFileId += "_MCSpring11";  
  else if( inFileListName.Contains("Summer11") ) outFileId += "_MCSummer11";  
  else outFileId += "_MC";  

  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    newFilesEtaPtAve[etaBin] = std::vector<TFile*>(binAdmin.nPtBins(hlt,etaBin));
    newTreesEtaPtAve[etaBin] = std::vector<TTree*>(binAdmin.nPtBins(hlt,etaBin));
    if( writeAllTrees ) {
      newFilesEtaPt[etaBin] = std::vector<TFile*>(binAdmin.nPtBins(hlt,etaBin));
      newTreesEtaPt[etaBin] = std::vector<TTree*>(binAdmin.nPtBins(hlt,etaBin));
      newFilesEtaGenPtGenAve[etaBin] = std::vector<TFile*>(binAdmin.nPtBins(hlt,etaBin));
      newTreesEtaGenPtGenAve[etaBin] = std::vector<TTree*>(binAdmin.nPtBins(hlt,etaBin));
    }
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
      TString etaPtId = "_Eta"+util::toTString(etaBin)+"_Pt"+util::toTString((binAdmin.hltMinPtBin(hlt,etaBin)+ptBin));
      
      TString name = outFilePrefix+"_PtAveBins"+outFileId+etaPtId+".root";
      newFilesEtaPtAve[etaBin][ptBin] = new TFile(name,"RECREATE");
      newTreesEtaPtAve[etaBin][ptBin] = oldChain->CloneTree(0);
      newTreesEtaPtAve[etaBin][ptBin]->SetDirectory(newFilesEtaPtAve[etaBin][ptBin]);
      
      if( writeAllTrees ) {
	name = outFilePrefix+"_PtBins"+outFileId+etaPtId+".root";
	newFilesEtaPt[etaBin][ptBin] = new TFile(name,"RECREATE");
	newTreesEtaPt[etaBin][ptBin] = oldChain->CloneTree(0);
	newTreesEtaPt[etaBin][ptBin]->SetDirectory(newFilesEtaPt[etaBin][ptBin]);
	
	name = outFilePrefix+"_PtGenAveBins"+outFileId+etaPtId+".root";
	newFilesEtaGenPtGenAve[etaBin][ptBin] = new TFile(name,"RECREATE");
	newTreesEtaGenPtGenAve[etaBin][ptBin] = oldChain->CloneTree(0);
	newTreesEtaGenPtGenAve[etaBin][ptBin]->SetDirectory(newFilesEtaGenPtGenAve[etaBin][ptBin]);
      }
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
      newTreesEtaPtAve[etaBin][ptBin]->Branch("L2L3CorrJetColJetIdx",corrJetIdx,"L2L3CorrJetColJetIdx[NobjJet]/I");
      newTreesEtaPtAve[etaBin][ptBin]->Branch("WeightSpectrum",&wSpec,"WeightSpectrum/F");
      if( writeAllTrees ) {
	newTreesEtaPt[etaBin][ptBin]->Branch("L2L3CorrJetColJetIdx",corrJetIdx,"L2L3CorrJetColJetIdx[NobjJet]/I");
	newTreesEtaPt[etaBin][ptBin]->Branch("WeightSpectrum",&wSpec,"WeightSpectrum/F");
	newTreesEtaGenPtGenAve[etaBin][ptBin]->Branch("L2L3CorrJetColJetIdx",corrJetIdx,"L2L3CorrJetColJetIdx[NobjJet]/I");
	newTreesEtaGenPtGenAve[etaBin][ptBin]->Branch("WeightSpectrum",&wSpec,"WeightSpectrum/F");
      }
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
       std::cerr << "WARNING: nObjJet = " << nObjJet << " > " << maxNJet << ". Skipping event!\n";
       ++nMaxNJet;
       continue;
     }


     if( isData ) {		// HLT selection
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
     } else {			// PU re-weighting

       wSpec = Weight;
       if( PUMCNumVtx < static_cast<int>(weightsPU.size()) ) {
 	wPU = weightsPU.at(PUMCNumVtx);
       } else {
 	std::cerr << "ERROR: no event weights for PUMCNumVtx = " << PUMCNumVtx << std::endl;
 	wPU = weightsPU.back();
       }

       // PU reweighting
       Weight *= wPU;
     }


     // Order corrected jets
     corrJets.clear();
     for(int j = 0; j < nObjJet; ++j) {

       //       // Vary JES
       //       double scale = 1. - JetCorrUncert[j];
       //       if( scale == scale ) {
       //  	//jetCorrL2L3[j] *= scale;
       //  	JetCorrL3[j] *= scale;
       //       } else {
       //  	std::cerr << "ERROR: JetCorrUncert[" << j << "] NAN" << std::endl;
       //       }

//        // Remove residual correction
//        jetCorrL2L3[j] = JetCorrL2[j]*JetCorrL3[j];

       corrJets.add(j,jetCorrL1[j]*jetCorrL2L3[j]*jetPt[j]);
     }
     corrJets.sort();

     // Set branch corrected jet indices
     for(int j = 0; j < nObjJet; ++j) {
       corrJetIdx[j] = corrJets(j);
     }

   

     // Data driven dijet selection
     if( nObjJet < 2 ) {
       ++nDijets;
     } else if( std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets(0)]-jetPhi[corrJets(1)])) < minDeltaPhi ) {
       ++nDeltaPhi;
     } else if( !jetID[corrJets(0)] || !jetID[corrJets(1)] ) {
       ++nJetID;
     } else {
       unsigned int etaBin = 1000;
       if( binAdmin.findSameEtaBin(jetEta[corrJets(0)],jetEta[corrJets(1)],etaBin) ) {
	
    	// Find ptAve bin
    	double ptAve = 0.5*(corrJets.pt(0)+corrJets.pt(1));
    	unsigned int ptAveBin = 1000;
    	if( binAdmin.findPtBin(hlt,ptAve,etaBin,ptAveBin) ) {
    	  ptAveBin -= binAdmin.hltMinPtBin(hlt,etaBin);
    	  newTreesEtaPtAve[etaBin][ptAveBin]->Fill();
    	} else {
    	  ++nPtAve;
    	}

  	if( writeAllTrees ) {
  	  // Find pt bin
  	  // Selection A:
  	  // - jet1: defines pt bin, i.e. ptMin < pt(jet1) < ptMax
  	  // - jet2: just minimum pt requirement, i.e. pt(jet2) > 45 GeV
  	  if( corrJets.pt(1) > 45. ) {
  	    unsigned int ptBin = 1000;
  	    if( binAdmin.findPtBin(hlt,corrJets.pt(0),etaBin,ptBin) ) {
  	      ptBin -= binAdmin.hltMinPtBin(hlt,etaBin);
  	      newTreesEtaPt[etaBin][ptBin]->Fill();
  	    }
  	  }
  	  // Selection B:
  	  // - jet2: defines pt bin, i.e. ptMin < pt(jet2) < ptMax
  	  // - jet1: just minimum pt requirement, i.e. pt(jet1) > 45 GeV
  	  if( corrJets.pt(0) > 45. ) {
  	    unsigned int ptBin = 1000;
  	    if( binAdmin.findPtBin(hlt,corrJets.pt(1),etaBin,ptBin) ) {
  	      ptBin -= binAdmin.hltMinPtBin(hlt,etaBin);
  	      corrJetIdx[0] = corrJets(1);
  	      corrJetIdx[1] = corrJets(0);
  	      newTreesEtaPt[etaBin][ptBin]->Fill();
  	      corrJetIdx[0] = corrJets(0);
  	      corrJetIdx[1] = corrJets(1);
  	    }
  	  }
  	}
       } else {
   	++nEta;
       }
     }
     

     if( !isData ) {
       
       // GenJet dijet selection
       if( NobjGenJet > 1 ) {
	 unsigned int etaGenBin = 1000;
	 if( binAdmin.findSameEtaBin(GenJetColEta[0],GenJetColEta[1],etaGenBin) ) {
	   // Indices of reco-jets matched to gen-jets: the reco-jet closest
	   // in deltaR is used; this matching has already been done
	   // when filling the original ntuples
	   int recoIdx0 = GenJetColJetIdx[0];
	   int recoIdx1 = GenJetColJetIdx[1];
	   // Kinematics
	   double ptGen0 = GenJetColPt[0];
	   double ptGen1 = GenJetColPt[1];
	   double ptGenAve = 0.5*(ptGen0+ptGen1);
	   double ptGen2 = NobjGenJet>2 ? GenJetColPt[2] : 0.;
	   // Response of leading two gen-jets
	   double r0 = 0.;
	   double r1 = 0.;
	   if( ptGen0 > 0. && ptGen1 > 0. ) {
	     r0 = jetCorrL1[recoIdx0]*jetCorrL2L3[recoIdx0]*jetPt[recoIdx0]/ptGen0;
	     r1 = jetCorrL1[recoIdx1]*jetCorrL2L3[recoIdx1]*jetPt[recoIdx1]/ptGen1;
	   }
	   // DeltaR between gen-jet and matched reco-jet
	   double dr0 = util::deltaR(jetEta[recoIdx0],GenJetColEta[0],jetPhi[recoIdx0],GenJetColPhi[0]);
	   double dr1 = util::deltaR(jetEta[recoIdx1],GenJetColEta[1],jetPhi[recoIdx1],GenJetColPhi[1]);
	   bool closeMatch = ( dr0 < maxDeltaR && dr1 < maxDeltaR );
	   // Do the matched reco-jets pass the jet id?
	   bool passJetID = ( jetID[recoIdx0] && jetID[recoIdx1] );
	   
	   // MC truth response only for good, closely matched reco-jets
	   if( passJetID ) {
	     hRespVsPtGenLogNoDeltaR.at(etaGenBin)->Fill(ptGen0,r0);
	     hRespVsPtGenLogNoDeltaR.at(etaGenBin)->Fill(ptGen1,r1);
	     if( closeMatch ) {
	       hRespVsPtGenLogNoDeltaPhi.at(etaGenBin)->Fill(ptGen0,r0);
	       hRespVsPtGenLogNoDeltaPhi.at(etaGenBin)->Fill(ptGen1,r1);
	     }
	   }
	   
	   // Require back-to-back gen-jets
	   if( std::abs(TVector2::Phi_mpi_pi(GenJetColPhi[0]-GenJetColPhi[1])) > minDeltaPhi ) {
	     // Fill various MC truth response distributions
	     // MC truth response only for good, closely matched reco-jets
	     if( passJetID && closeMatch ) {
	       hRespVsPtGen.at(etaGenBin)->Fill(ptGen0,r0);
	       hRespVsPtGen.at(etaGenBin)->Fill(ptGen1,r1);
	       hRespVsPtGenLog.at(etaGenBin)->Fill(ptGen0,r0);
	       hRespVsPtGenLog.at(etaGenBin)->Fill(ptGen1,r1);
	       hRespVsPtGenPUReweighted.at(etaGenBin)->Fill(ptGen0,r0,wPU);
	       hRespVsPtGenPUReweighted.at(etaGenBin)->Fill(ptGen1,r1,wPU);
	       hRespVsPtGenLogPUReweighted.at(etaGenBin)->Fill(ptGen0,r0,wPU);
	       hRespVsPtGenLogPUReweighted.at(etaGenBin)->Fill(ptGen1,r1,wPU);
	       if( PUMCNumVtx < 5 ) {
		 hRespVsPtGenLogPULow.at(etaGenBin)->Fill(ptGen0,r0);
		 hRespVsPtGenLogPULow.at(etaGenBin)->Fill(ptGen1,r1);
	       } else if( PUMCNumVtx > 9 ) {
		 hRespVsPtGenLogPUHigh.at(etaGenBin)->Fill(ptGen0,r0);
		 hRespVsPtGenLogPUHigh.at(etaGenBin)->Fill(ptGen1,r1);
	       }
	     }
	      
	     // Pt soft cuts
	     for(int ptSoftBin = binAdmin.nPtSoftBins()-1; ptSoftBin >= 0; --ptSoftBin) {
	       if( ptGen2 > binAdmin.ptSoftMax(ptSoftBin)*ptGenAve ) break;
	       // Gen-jet spectrum independent from reco-jet properties
	       hPtGen.at(etaGenBin).at(ptSoftBin)->Fill(ptGen0,wSpec);
	       hPtGen.at(etaGenBin).at(ptSoftBin)->Fill(ptGen1,wSpec);
	       // MC truth response only for good, closely matched reco-jets
	       if( passJetID && closeMatch ) {
		 hRespVsPtGenLogVsPtSoftPUReweighted.at(etaGenBin).at(ptSoftBin)->Fill(ptGen0,r0,wSpec*wPU);
		 hRespVsPtGenLogVsPtSoftPUReweighted.at(etaGenBin).at(ptSoftBin)->Fill(ptGen1,r1,wSpec*wPU);
	       }
	     }
	      
	     if( writeAllTrees ) {
	       // Find ptGenAve bin
	       unsigned int ptGenAveBin = 1000;
	       if( binAdmin.findPtBin(ptGenAve,etaGenBin,ptGenAveBin) ) {
		 for(int j = 0; j < nObjJet; ++j) {
		   corrJetIdx[j] = GenJetColJetIdx[j];
		 }
		 newTreesEtaGenPtGenAve[etaGenBin][ptGenAveBin]->Fill();
	       }
	     }
	   }
	 }
       }
     }
   
     
     if( i%50000 == 0 ) {
       for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
  	for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
 	  newTreesEtaPtAve[etaBin][ptBin]->AutoSave();
 	  if( writeAllTrees ) {
 	    newTreesEtaPt[etaBin][ptBin]->AutoSave();
 	    newTreesEtaGenPtGenAve[etaBin][ptBin]->AutoSave();
 	  }
  	}
       }
     }
   }     // End of loop over entries

     for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
       for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
 	newTreesEtaPtAve[etaBin][ptBin]->AutoSave();
 	if( writeAllTrees ) {
 	  newTreesEtaPt[etaBin][ptBin]->AutoSave();
 	  newTreesEtaGenPtGenAve[etaBin][ptBin]->AutoSave();
 	}
       }
     }

     if( !isData ) {

       TFile fileResponse(outFilePath+"/Kalibri_MCTruthResponse.root","RECREATE");
       TFile fileSpetrum(outFilePath+"/Kalibri_MCTruthSpectrum.root","RECREATE");
       for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
 	fileResponse.WriteTObject(hRespVsPtGen.at(etaBin));
 	fileResponse.WriteTObject(hRespVsPtGenLog.at(etaBin));
	fileResponse.WriteTObject(hRespVsPtGenLogNoDeltaR.at(etaBin));
 	fileResponse.WriteTObject(hRespVsPtGenLogNoDeltaPhi.at(etaBin));
 	fileResponse.WriteTObject(hRespVsPtGenPUReweighted.at(etaBin));
 	fileResponse.WriteTObject(hRespVsPtGenLogPUReweighted.at(etaBin));
 	fileResponse.WriteTObject(hRespVsPtGenLogPULow.at(etaBin));
 	fileResponse.WriteTObject(hRespVsPtGenLogPUHigh.at(etaBin));

 	for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
 	  fileResponse.WriteTObject(hRespVsPtGenLogVsPtSoftPUReweighted.at(etaBin).at(ptSoftBin));
 	  fileSpetrum.WriteTObject(hPtGen.at(etaBin).at(ptSoftBin));
 	}
       }
       fileResponse.Close();
       fileSpetrum.Close();
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
 	std::cout << "  " << newTreesEtaPtAve[etaBin][ptBin]->GetEntries() << " events with " << binAdmin.etaMin(etaBin) << " < |eta(1,2)| < " << binAdmin.etaMax(etaBin) << " to file '" << newFilesEtaPtAve[etaBin][ptBin]->GetName() << "'" << std::endl;
       }
     }



     // ++++ Clean up +++++++++++++++++++++++++++++++++++++++++++++

     for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
       for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
 	newFilesEtaPtAve[etaBin][ptBin]->Close();
 	delete newFilesEtaPtAve[etaBin][ptBin];

 	if( writeAllTrees ) {
 	  newFilesEtaPt[etaBin][ptBin]->Close();
 	  newFilesEtaGenPtGenAve[etaBin][ptBin]->Close();
 	  delete newFilesEtaPt[etaBin][ptBin];
 	  delete newFilesEtaGenPtGenAve[etaBin][ptBin];
 	}
       }
     }
     delete oldChain;

     if( !isData ) {
       for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
 	delete hRespVsPtGen.at(etaBin);
 	delete hRespVsPtGenLog.at(etaBin);
 	delete hRespVsPtGenLogNoDeltaR.at(etaBin);
 	delete hRespVsPtGenLogNoDeltaPhi.at(etaBin);
 	delete hRespVsPtGenPUReweighted.at(etaBin);
 	delete hRespVsPtGenLogPUReweighted.at(etaBin);
 	delete hRespVsPtGenLogPULow.at(etaBin);
 	delete hRespVsPtGenLogPUHigh.at(etaBin);
  
 	for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
 	  delete hRespVsPtGenLogVsPtSoftPUReweighted.at(etaBin).at(ptSoftBin);
 	  delete hPtGen.at(etaBin).at(ptSoftBin);
 	}
       }
     }
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
