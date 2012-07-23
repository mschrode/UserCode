// $Id: writeDijetSkims.C,v 1.11 2012/02/04 21:49:09 mschrode Exp $
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


enum MCSampleType { Flat, PtHatBinned };

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


// Generate weights for different PU scenarios for given
// data PU distribution
// Code from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
// S3 = Flat10 (Summer11)
// S6 = Flat25? (Fall11)
// --------------------------------------------------
std::vector<double> generate_PUS3_weights(const TH1* data_npu_estimated) {
  std::cout << "Generating PU weights for S3 scenario (Summer11)" << std::endl;
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

std::vector<double> generate_PUS6_weights(const TH1* data_npu_estimated) {
  std::cout << "Generating PU weights for S6 scenario (Fall11)" << std::endl;
  const double npu_probs[50] = {
    0.003388501,
    0.010357558,
    0.024724258,
    0.042348605,
    0.058279812,
    0.068851751,
    0.072914824,
    0.071579609,
    0.066811668,
    0.060672356,
    0.054528356,
    0.04919354,
    0.044886042,
    0.041341896,
    0.0384679,
    0.035871463,
    0.03341952,
    0.030915649,
    0.028395374,
    0.025798107,
    0.023237445,
    0.020602754,
    0.0180688,
    0.015559693,
    0.013211063,
    0.010964293,
    0.008920993,
    0.007080504,
    0.005499239,
    0.004187022,
    0.003096474,
    0.002237361,
    0.001566428,
    0.001074149,
    0.000721755,
    0.000470838,
    0.00030268,
    0.000184665,
    0.000112883,
    6.74043E-05,
    3.82178E-05,
    2.22847E-05,
    1.20933E-05,
    6.96173E-06,
    3.4689E-06,
    1.96172E-06,
    8.49283E-07,
    5.02393E-07,
    2.15311E-07,
    9.56938E-08
  };
  std::vector<double> result(50);
  double s = 0.0;
  for(int npu = 0; npu < 50; ++npu) {
    double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
    result[npu] = npu_estimated / npu_probs[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(int npu = 0; npu < 50; ++npu) {
    result[npu] /= s;
  }
  return result;
}


// Get weight for ptHat binned MC samples
// The weight is set for an integrated luminosity
// of 100 / pb
//
// The weights returned are hard-coded for the
// Summer11 PYTHIA QCD samples Tune Z2,
// /QCD_Pt-*to*_TuneZ2_7TeV_pythia6/Summer11-PU_S3_START42_V11-v2/AODSIM
// --------------------------------------------------
double getWeightForPtHatBinnedMCSample(double ptHat) {
  double lumi = 100.;
  unsigned int numEvts = 0;
  double xs = 0.;
  double scale = 1.;		// To take into account missing jobs
  if( ptHat < 15. ) {
    numEvts   = 1650000;
    xs        = 3.675e+10;
  } else if( ptHat < 30. ) {
    numEvts   = 11000000;
    xs        = 8.159e+08;
    scale = 382./381.;
  } else if( ptHat < 50. ) {
    numEvts   = 6583068;
    xs        = 5.312e+07;
  } else if( ptHat < 80. ) {
    numEvts   = 6600000;
    xs        = 6.359e+06;
    scale = 233./232.;
  } else if( ptHat < 120. ) {
    numEvts   = 6589956;
    xs        = 7.843e+05;
  } else if( ptHat < 170. ) {
    numEvts   = 6127528;
    xs        = 1.151e+05;
  } else if( ptHat < 300. ) {
    numEvts   = 6220160;
    xs        = 2.426e+04;
    scale = 210./209.;
  } else if( ptHat < 470. ) {
    numEvts   = 6432669;
    xs        = 1.168e+03;
    scale = 217./215.;
  } else if( ptHat < 600. ) {
    numEvts   = 3990085;
    xs        = 7.022e+01;
    scale     = 134./129.;
  } else if( ptHat < 800. ) {
    numEvts   = 4245695;
    xs        = 1.555e+01;
  } else if( ptHat < 1000. ) {
    numEvts   = 4053888;
    xs        = 1.844e+00;
    scale     = 137./134.;
  } else if( ptHat < 1400. ) {
    numEvts   = 2093222;
    xs        = 3.321e-01;
  } else if( ptHat < 1800. ) {
    numEvts   = 2196200;
    xs        = 1.087e-02;
    scale     = 74./73.;
  } else {
    numEvts   = 293139;
    xs        = 3.575e-04;
  }

  return numEvts > 0 ? scale * lumi * xs / numEvts : 0.;
}



// --------------------------------------------------
void writeDijetSkims(bool isData, unsigned int maxHltThres = 0) {
  
  
  // ++++ Set parameters +++++++++++++++++++++++++++++++++++++++
  
  std::cout << "Setting parameters" << std::endl;
  
  //  const TString inFileListName = "input/Analysis2011/Kalibri_JetRun2011A_V3_163337-167151_L1FastJet_AK5PF";
  //const TString inFileListName = "input/Analysis2011/Kalibri_MCSummer11_QCDFlat_PythiaZ2_PUS3_L1FastJet_AK5PF";
  //  const TString inFileListName = "input/Analysis2011/Kalibri_JetRun2011_AK5PF_L1FastJet_V10";
  //  const TString inFileListName = "input/Analysis2011/Kalibri_MCFall11_QCDFlat_PythiaZ2_PUS6_AK5PF_L1FastJet_V10";
  const TString inFileListName = "input/Analysis2011/Kalibri_WJetsToLNu_TuneZ2_7TeV-madgraph-tauola-Fall11-PU_S6_START42_AK5PF_L1FastJet_V10";

  const TString outFilePath = "~/lustre/tmp/test";
  const TString configAdm = "../resolutionFit/config/Analysis2011/Binning/BinningAdmin2011_v2.cfg";
  const TString configBin = "../resolutionFit/config/Analysis2011/Binning/Binning2011_v2_skims.cfg";

  const MCSampleType mcSampleType_ = Flat;
  const bool writeAllTrees = true;
  
  const int nEvts = -10000000;
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
  unsigned int nJetID = 0;
  unsigned int nLargeResp = 0;
  unsigned int nDeltaPhi = 0;
  unsigned int nEta = 0;
  unsigned int nPtAve = 0;
  
  // Event weights
  double wPU = 1.;		// PU weight
  float wSpec = 1.;		// Spectrum weight (for flat samples, taken from input ntuple)
  
  

  // ++++ Prepare objects ++++++++++++++++++++++++++++++++++++++

  std::cout << "Defining binning" << std::endl;

  const sampleTools::BinningAdmin binAdmin(configAdm,configBin);
  const TString hlt = isData ? binAdmin.triggerName(maxHltThres) : "none";
  binAdmin.printBinning();
  if( isData ) binAdmin.print(hlt);

  std::vector<double> weightsPU;
  // MC-truth response
  std::vector<TH2*> hRespVsPtGen; // Default selection
  std::vector<TH2*> hRespVsPtGenUncorrected; // Default selection
  std::vector<TH2*> hRespVsPtGenUncorrectedLowPU; // Default selection
  std::vector<TH2*> hRespVsPtGenL1Corrected; // Default selection
  TH1* hPtGenBinsForRespVsEtaGen = 0;
  std::vector<TH2*> hRespVsEtaGen; // Default selection
  std::vector<TH2*> hRespVsEtaGenUncorrected; // Default selection
  std::vector<TH2*> hRespVsEtaGenL1Corrected; // Default selection
  std::vector<TH2*> hRespVsPtGenPUMC; // MC pile-up scenario
  std::vector<TH2*> hRespVsPtGenPULess05;
  std::vector<TH2*> hRespVsPtGenPULess10;
  std::vector<TH2*> hRespVsPtGenPULess15;
  std::vector<TH2*> hRespVsPtGenPULess99;
  std::vector<TH2*> hRespVsPtGenNoJetID;
  std::vector<TH2*> hRespVsPtGenNoDeltaR;
  std::vector<TH2*> hRespVsPtGenDeltaRLess05;
  std::vector<TH2*> hRespVsPtGenDeltaRLess10;
  std::vector<TH2*> hRespVsPtGenDeltaRLess15;
  std::vector<TH2*> hRespVsPtGenDeltaRLess20;
  std::vector<TH2*> hRespVsPtGenDeltaRLess25;
  std::vector<TH2*> hRespVsPtGenDeltaRLess30;
  std::vector<TH2*> hDeltaRVsPtGen;
  std::vector< std::vector<TH2*> > hRespVsPtGenVsPtSoft; // [etaBin][pt3Bin]
  std::vector< std::vector<TH2*> > hRespVsPtGenVsPdgId; // [etaBin][pt3Bin]
  std::vector< std::vector<TH2*> > hRespVsPtGenUncorrectedVsPdgId; // [etaBin][pt3Bin]
  std::vector< std::vector<TH2*> > hRespVsPtGenVsPdgIdDijets; // [etaBin][pt3Bin]
  // PtGen spectrum
  std::vector< std::vector<TH1*> > hPtGen; // [etaBin][pt3Bin]
  std::vector< std::vector<TH1*> > hPtGenUnweighted; // [etaBin][pt3Bin]

  // N-1 and control plots
  std::vector< std::vector<TH1*> > hEta;
  std::vector< std::vector<TH1*> > hEtaNMin1;
  std::vector< std::vector<TH1*> > hDeltaPhiNMin1;
  std::vector< std::vector<TH1*> > hPtRelNMin1;
  std::vector< std::vector<TH2*> > hDeltaPhiVsPt3Rel;
  std::vector< std::vector<TH2*> > hPt3RelGenVsReco;
  std::vector< std::vector<TH2*> > hPt3RelRecoVsImbalGen;
  std::vector< std::vector<TH2*> > hPt3RelGenVsImbalGen;
  std::vector< std::vector<TH1*> > hPtJet;
  std::vector<TH1*> hPtAve;
  
  if( !isData ) {
    std::cout << "Preparing PU reweighting" << std::endl;
    
//     TH1* hDataPU = util::FileOps::readTH1("~/PileUp/Pileup_2011_to_172255.root","pileup");
//     weightsPU = generate_PUS3_weights(hDataPU);
    TH1* hDataPU = util::FileOps::readTH1("~/PileUp/ExpectedPileUpDist_160404_180252_SelfCombined.root","pileup");
    weightsPU = generate_PUS6_weights(hDataPU);
    
    std::cout << "Preparing histograms" << std::endl;
    
    // MC-truth histograms

    // Equidistant log bins
    std::vector<double> binEdges(26);
    util::HistOps::equidistLogBins(binEdges,binEdges.size()-1,4.,2700.);

    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      TString title = util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin));
      
      TH2* h = new TH2D("hRespVsPtGen_Eta"+util::toTString(etaBin),"p^{gen}_{T} (GeV);Response",
			binEdges.size()-1,&(binEdges.front()),201,0.,2.);
      h->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, PU reweighted");
      h->SetNdivisions(505);
      h->Sumw2();
      hRespVsPtGen.push_back(h);

      hRespVsPtGenUncorrected.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenUncorrected_Eta"+util::toTString(etaBin))));
      hRespVsPtGenUncorrected.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, Uncorrected");

      hRespVsPtGenUncorrectedLowPU.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenUncorrectedLowPU_Eta"+util::toTString(etaBin))));
      hRespVsPtGenUncorrectedLowPU.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, Uncorrected, N(PU)<2");

      hRespVsPtGenL1Corrected.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenL1Corrected_Eta"+util::toTString(etaBin))));
      hRespVsPtGenL1Corrected.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, Uncorrected, L1 corrected");

      hRespVsPtGenPUMC.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenPUMC_Eta"+util::toTString(etaBin))));
      hRespVsPtGenPUMC.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, PU MC");

      hRespVsPtGenPULess05.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenPULess05_Eta"+util::toTString(etaBin))));
      hRespVsPtGenPULess05.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, N(PU)<5");

      hRespVsPtGenPULess10.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenPULess10_Eta"+util::toTString(etaBin))));
      hRespVsPtGenPULess10.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, 4<N(PU)<10");

      hRespVsPtGenPULess15.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenPULess15_Eta"+util::toTString(etaBin))));
      hRespVsPtGenPULess15.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, 9<N(PU)<15");

      hRespVsPtGenPULess99.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenPULess99_Eta"+util::toTString(etaBin))));
      hRespVsPtGenPULess99.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, N(PU)>14");

      hRespVsPtGenNoJetID.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenNoJetID_Eta"+util::toTString(etaBin))));
      hRespVsPtGenNoJetID.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", PU reweighted");

      hRespVsPtGenNoDeltaR.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenNoDeltaR_Eta"+util::toTString(etaBin))));
      hRespVsPtGenNoDeltaR.back()->SetTitle(title+"JetID, PU reweighted");

      hRespVsPtGenDeltaRLess05.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenDeltaRLess05_Eta"+util::toTString(etaBin))));
      hRespVsPtGenDeltaRLess05.back()->SetTitle(title+", |#DeltaR|<0.05, JetID, PU reweighted");

      hRespVsPtGenDeltaRLess10.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenDeltaRLess10_Eta"+util::toTString(etaBin))));
      hRespVsPtGenDeltaRLess10.back()->SetTitle(title+", |#DeltaR|<0.05, JetID, PU reweighted");

      hRespVsPtGenDeltaRLess15.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenDeltaRLess15_Eta"+util::toTString(etaBin))));
      hRespVsPtGenDeltaRLess15.back()->SetTitle(title+", |#DeltaR|<0.15, JetID, PU reweighted");

      hRespVsPtGenDeltaRLess20.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenDeltaRLess20_Eta"+util::toTString(etaBin))));
      hRespVsPtGenDeltaRLess20.back()->SetTitle(title+", |#DeltaR|<0.20, JetID, PU reweighted");

      hRespVsPtGenDeltaRLess25.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenDeltaRLess25_Eta"+util::toTString(etaBin))));
      hRespVsPtGenDeltaRLess25.back()->SetTitle(title+", |#DeltaR|<0.25, JetID, PU reweighted");

      hRespVsPtGenDeltaRLess30.push_back(static_cast<TH2*>(h->Clone("hRespVsPtGenDeltaRLess30_Eta"+util::toTString(etaBin))));
      hRespVsPtGenDeltaRLess30.back()->SetTitle(title+", |#DeltaR|<0.30, JetID, PU reweighted");

      h = new TH2D("hDeltaRVsPtGen_Eta"+util::toTString(etaBin),"p^{gen}_{T} (GeV);#DeltaR",
		   binEdges.size()-1,&(binEdges.front()),50,0.,1.);
      h->SetTitle(title+", JetID");
      h->SetNdivisions(505);
      h->Sumw2();
      hDeltaRVsPtGen.push_back(h);

      std::vector<TH2*> vtmp;
      for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
	title = util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin))+", p^{rel}_{T,3} < "+util::toTString(binAdmin.ptSoftMax(ptSoftBin))+", |#DeltaR|<"+util::toTString(maxDeltaR)+", |#Delta#Phi|>"+util::toTString(minDeltaPhi)+", PU reweighted;p^{gen}_{T} (GeV);Response";
	h = new TH2D("hRespVsPtGenVsPtSoft_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin),title+";p^{gen}_{T} (GeV);Response",binEdges.size()-1,&(binEdges.front()),201,0.,2.);
	h->SetNdivisions(505);
	h->Sumw2();
	vtmp.push_back(h);
      }
      hRespVsPtGenVsPtSoft.push_back(vtmp);
      
      vtmp.clear();
      for(unsigned int id = 0; id < 22; ++id) {
	title = util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin))+", PDG Id = "+util::toTString(id)+", |#DeltaR|<"+util::toTString(maxDeltaR)+", PU reweighted;p^{gen}_{T} (GeV);Response";
	h = new TH2D("hRespVsPtGenVsPdgId_Eta"+util::toTString(etaBin)+"_PdgId"+util::toTString(id),title+";p^{gen}_{T} (GeV);Response",binEdges.size()-1,&(binEdges.front()),201,0.,2.);
	h->SetNdivisions(505);
	h->Sumw2();
	vtmp.push_back(h);
      }
      hRespVsPtGenVsPdgId.push_back(vtmp);

      vtmp.clear();
      for(unsigned int id = 0; id < 22; ++id) {
	title = util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin))+", PDG Id = "+util::toTString(id)+", |#DeltaR|<"+util::toTString(maxDeltaR)+", |#Delta#Phi| > 2.7, p^{rel}_{T,3} < 0.14"+", PU reweighted, Uncorrected;p^{gen}_{T} (GeV);Response";
	h = new TH2D("hRespVsPtGenUncorrectedVsPdgId_Eta"+util::toTString(etaBin)+"_PdgId"+util::toTString(id),title+";p^{gen}_{T} (GeV);Response",binEdges.size()-1,&(binEdges.front()),201,0.,2.);
	h->SetNdivisions(505);
	h->Sumw2();
	vtmp.push_back(h);
      }
      hRespVsPtGenUncorrectedVsPdgId.push_back(vtmp);

      vtmp.clear();
      for(unsigned int id = 0; id < 22; ++id) {
	title = util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin))+", PDG Id = "+util::toTString(id)+", |#DeltaR|<"+util::toTString(maxDeltaR)+", |#Delta#Phi| > 2.7, p^{rel}_{T,3} < 0.14"+", PU reweighted;p^{gen}_{T} (GeV);Response";
	h = new TH2D("hRespVsPtGenVsPdgIdDijets_Eta"+util::toTString(etaBin)+"_PdgId"+util::toTString(id),title+";p^{gen}_{T} (GeV);Response",binEdges.size()-1,&(binEdges.front()),201,0.,2.);
	h->SetNdivisions(505);
	h->Sumw2();
	vtmp.push_back(h);
      }
      hRespVsPtGenVsPdgIdDijets.push_back(vtmp);
    }

    binEdges.clear();
    binEdges = std::vector<double>(11);
    util::HistOps::equidistLogBins(binEdges,binEdges.size()-1,4.,2700.);
    hPtGenBinsForRespVsEtaGen = new TH1D("hPtGenBinsForRespVsEtaGen",";p^{gen}_{T} (GeV)",binEdges.size()-1,&(binEdges.front()));
    for(int ptGenBin = 0; ptGenBin < hPtGenBinsForRespVsEtaGen->GetNbinsX(); ++ptGenBin) {
      double ptMin = hPtGenBinsForRespVsEtaGen->GetXaxis()->GetBinLowEdge(1+ptGenBin);
      double ptMax = hPtGenBinsForRespVsEtaGen->GetXaxis()->GetBinUpEdge(1+ptGenBin);
      TString title = util::toTString(ptMin,1)+" < p^{gen}_{T} < "+util::toTString(ptMax,1)+" GeV";
      
      TH2* h = new TH2D("hRespVsEtaGen_PtGen"+util::toTString(ptGenBin),";#eta^{gen};Response",52,-5.2,5.2,201,0.,2.);
      h->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, PU reweighted");
      h->SetNdivisions(505);
      h->Sumw2();
      hRespVsEtaGen.push_back(h);

      hRespVsEtaGenUncorrected.push_back(static_cast<TH2*>(h->Clone("hRespVsEtaGenUncorrected_PtGen"+util::toTString(ptGenBin))));
      hRespVsEtaGenUncorrected.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, Uncorrected");

      hRespVsEtaGenL1Corrected.push_back(static_cast<TH2*>(h->Clone("hRespVsEtaGenL1Corrected_PtGen"+util::toTString(ptGenBin))));
      hRespVsEtaGenL1Corrected.back()->SetTitle(title+", |#DeltaR|<"+util::toTString(maxDeltaR)+", JetID, L1 Corrected");
    }


    // PtGen spectra
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      std::vector<TH1*> vtmp1;
      std::vector<TH1*> vtmp2;
      for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
	TString title = util::toTString(binAdmin.etaMin(etaBin))+" < |#eta| < "+util::toTString(binAdmin.etaMax(etaBin))+", p^{rel}_{T,3} < "+util::toTString(binAdmin.ptSoftMax(ptSoftBin))+";p^{gen}_{T} (GeV);Events";

	TH1* h = new TH1D("hPtGen_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin),
			  title,500,0.,2000.);
	h->GetXaxis()->SetNdivisions(505);
	h->SetLineWidth(2);
	h->Sumw2();
	vtmp1.push_back(h);

	vtmp2.push_back(static_cast<TH1*>(h->Clone("hPtGenUnweighted_Eta"+util::toTString(etaBin)+"_PtSoft"+util::toTString(ptSoftBin))));
      }
      hPtGen.push_back(vtmp1);
      hPtGenUnweighted.push_back(vtmp2);
    }
  }

  // N-1 and control plots

  // Loop over eta bins
  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    hPtAve.push_back(new TH1D("hPtAve_EtaBin"+util::toTString(etaBin),";p^{ave}_{T} (GeV)",1500,0.,1500.));
    hPtAve.back()->Sumw2();

    // Vs jet rank
    std::vector<TH1*> hPtJetTmp;
    for(unsigned int i = 0; i < 4; ++i) {
      hPtJetTmp.push_back(new TH1D("hPtJet"+util::toTString(i)+"_EtaBin"+util::toTString(etaBin),";p_{T,"+util::toTString(i)+"} (GeV)",1500,0.,1500.));
      hPtJetTmp.back()->Sumw2();
    } 
    hPtJet.push_back(hPtJetTmp);
    // End vs jet rank

    // Vs ptAve bin
    std::vector<TH1*> hEtaTmp;
    std::vector<TH1*> hEtaNMin1Tmp;
    std::vector<TH1*> hDeltaPhiNMin1Tmp;
    std::vector<TH1*> hPtRelNMin1Tmp;
    std::vector<TH2*> hDeltaPhiVsPt3RelTmp;
    std::vector<TH2*> hPt3RelGenVsRecoTmp;
    std::vector<TH2*> hPt3RelRecoVsImbalGenTmp;
    std::vector<TH2*> hPt3RelGenVsImbalGenTmp;

    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
      TString binId = "_EtaBin"+util::toTString(etaBin)+"_PtBin"+util::toTString(ptBin);
      hEtaTmp.push_back(new TH1D("hEta"+binId,"N-1;#eta",102,-5.1,5.1));
      hEtaTmp.back()->Sumw2();
      hEtaNMin1Tmp.push_back(static_cast<TH1*>(hEtaTmp.back()->Clone("hEtaNMin1"+binId)));
      hDeltaPhiNMin1Tmp.push_back(new TH1D("hDeltaPhiNMin1"+binId,"N-1;#Delta#phi",320,0.,3.2));
      hDeltaPhiNMin1Tmp.back()->Sumw2();
      hPtRelNMin1Tmp.push_back(new TH1D("hPtRelNMin1"+binId,"N-1;p^{rel}_{T,3}",120,0.,1.2));
      hPtRelNMin1Tmp.back()->Sumw2();
      hDeltaPhiVsPt3RelTmp.push_back(new TH2D("hDeltaPhiVsPt3Rel"+binId,";p^{rel}_{T,3};#Delta#phi",120,0.,1.2,400,0.,4.));
      hDeltaPhiVsPt3RelTmp.back()->Sumw2();
      hPt3RelGenVsRecoTmp.push_back(new TH2D("hPt3RelGenVsReco"+binId,";p^{rel}_{T,3};p^{gen,rel}_{T,3}",120,0.,1.2,120,0.,1.2));
      hPt3RelGenVsRecoTmp.back()->Sumw2();
      hPt3RelRecoVsImbalGenTmp.push_back(new TH2D("hPt3RelRecoVsImbalGen"+binId,";#alpha_{imbal};p^{rel}_{T,3}",120,0.,1.2,120,0.,1.2));
      hPt3RelRecoVsImbalGenTmp.back()->Sumw2();
      hPt3RelGenVsImbalGenTmp.push_back(new TH2D("hPt3RelGenVsImbalGen"+binId,";#alpha_{imbal};p^{gen,rel}_{T,3}",120,0.,1.2,120,0.,1.2));
      hPt3RelGenVsImbalGenTmp.back()->Sumw2();
    }
    hEta.push_back(hEtaTmp);
    hEtaNMin1.push_back(hEtaNMin1Tmp);
    hDeltaPhiNMin1.push_back(hDeltaPhiNMin1Tmp);
    hPtRelNMin1.push_back(hPtRelNMin1Tmp);
    hDeltaPhiVsPt3Rel.push_back(hDeltaPhiVsPt3RelTmp);
    hPt3RelGenVsReco.push_back(hPt3RelGenVsRecoTmp);
    hPt3RelRecoVsImbalGen.push_back(hPt3RelRecoVsImbalGenTmp);
    hPt3RelGenVsImbalGen.push_back(hPt3RelGenVsImbalGenTmp);
    // End vs ptAve bin
  }// End loop over eta bins


  std::cout << "Preparing trees" << std::endl;

  // Open Kalibri ntuples
  TChain *oldChain = createTChain(inFileListName);
  
  // Deactivate branches not needed
  oldChain->SetBranchStatus("Track*",0);
  oldChain->SetBranchStatus("NobjTrack",0);
  oldChain->SetBranchStatus("Tow*",0);
  oldChain->SetBranchStatus("GenPart*",0);
  oldChain->SetBranchStatus("GenPartId*",1);
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
  Float_t         GenEvtScale;
  Int_t           VtxN = 0;
  Int_t           PUMCNumVtx = 0;
  Float_t         Weight = 1.;
  Int_t           GenPartId_algo[maxNJet];   //[NobjJet]
  Int_t           GenPartId_phys[maxNJet];   //[NobjJet]


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
  oldChain->SetBranchAddress("GenEvtScale",&GenEvtScale);
  oldChain->SetBranchAddress("VtxN",&VtxN);
  oldChain->SetBranchAddress("PUMCNumVtx",&PUMCNumVtx);
  oldChain->SetBranchAddress("Weight",&Weight);
  oldChain->SetBranchAddress("GenPartId_algo",GenPartId_algo);
  oldChain->SetBranchAddress("GenPartId_phys",GenPartId_phys);


  // Create a new file + a clone of old tree in new file
  // 1) per eta and ptAve bin;
  // 2) per etaGen and ptGenAve bin.

  // Prepare files and tress
  std::vector< std::vector<TFile*> > newFilesEtaPtAve(binAdmin.nEtaBins());
  std::vector< std::vector<TTree*> > newTreesEtaPtAve(binAdmin.nEtaBins());
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
  else if( inFileListName.Contains("Fall11") ) outFileId += "_MCFall11";  
  else outFileId += "_MC";  

  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    newFilesEtaPtAve[etaBin] = std::vector<TFile*>(binAdmin.nPtBins(hlt,etaBin));
    newTreesEtaPtAve[etaBin] = std::vector<TTree*>(binAdmin.nPtBins(hlt,etaBin));
    if( writeAllTrees ) {
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


    if( isData ) {
      // ---- for data ------------------------------------
      // HLT selection
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
    } else {	
      // ---- for MC ----------------------------------------------
      // Spectrum re-weighting
      if( mcSampleType_ == Flat ) {
	// Correct weight already in ntuples
	//	wSpec = Weight;
	wSpec = 1.;
      } else if( mcSampleType_ == PtHatBinned ) {
	// Compute weight from ptHat information
	wSpec = getWeightForPtHatBinnedMCSample(GenEvtScale);
      } else {
	std::cout << "WARNING: No spectrum re-weighting applied!" << std::endl;
	wSpec = 1.;
      }

      // PU re-weighting
      if( PUMCNumVtx < static_cast<int>(weightsPU.size()) ) {
 	wPU = weightsPU.at(PUMCNumVtx);
      } else {
	std::cerr << "ERROR: no event weights for PUMCNumVtx = " << PUMCNumVtx << " - using weights for PUMCNumVtx = " << weightsPU.size()-1 << " instead" <<  std::endl;
 	wPU = weightsPU.back();
      }

      // Combined event weight: spectrum weight x PU weight
      Weight = wSpec * wPU;
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
    } else if( !jetID[corrJets(0)] || !jetID[corrJets(1)] ) {
      ++nJetID;
    } else if( !isData &&
	       ( GenJetPt[corrJets(0)] > 0. && GenJetPt[corrJets(1)] > 0. ) &&
	       ( corrJets.pt(0)/GenJetPt[corrJets(0)] > 2.5 || corrJets.pt(1)/GenJetPt[corrJets(1)] > 2.5 ) ) {
      ++nLargeResp;
      //std::cout << "WARNING: Large response, omitting event" << std::endl;
    } else {
      double deltaPhi = std::abs(TVector2::Phi_mpi_pi(jetPhi[corrJets(0)]-jetPhi[corrJets(1)]));
      bool passDeltaPhi = deltaPhi > minDeltaPhi;
      double weight = wPU*wSpec;
      double ptAve = 0.5*(corrJets.pt(0)+corrJets.pt(1));
      double pt2 = nObjJet > 2 ? corrJets.pt(2) : 0.;
      double ptrel = pt2/ptAve;
      bool passSoftJet = ptrel < 0.14;
      double eta0 = jetEta[corrJets(0)];
      double eta1 = jetEta[corrJets(1)];

      // N-1 eta selection
      unsigned int etaBin = 1000;
      unsigned int ptAveBin = 1000;
      if( passDeltaPhi && passSoftJet ) {
	if( binAdmin.findEtaBin(eta0,etaBin) ) {
	  if( binAdmin.findPtBin(hlt,ptAve,etaBin,ptAveBin) ) {
	    hEtaNMin1.at(etaBin).at(ptAveBin)->Fill(eta0,weight);
	  }
	}
	if( binAdmin.findEtaBin(eta1,etaBin) ) {
	  if( binAdmin.findPtBin(hlt,ptAve,etaBin,ptAveBin) ) {
	    hEtaNMin1.at(etaBin).at(ptAveBin)->Fill(eta1,weight);
	  }
	}
      }
       
      if( binAdmin.findSameEtaBin(eta0,eta1,etaBin) ) { // Same eta bin
	if( binAdmin.findPtBin(hlt,ptAve,etaBin,ptAveBin) ) { // PtAve bin
	   
	  hDeltaPhiVsPt3Rel.at(etaBin).at(ptAveBin)->Fill(ptrel,deltaPhi,weight);
	   
	  if( passSoftJet ) {
	    hDeltaPhiNMin1.at(etaBin).at(ptAveBin)->Fill(deltaPhi,weight);
	  }
	  if( passDeltaPhi ) {
	    hPtRelNMin1.at(etaBin).at(ptAveBin)->Fill(ptrel,weight);
	    if( !isData ) {
	      if( nObjJet > 2 ) {
		hPt3RelGenVsReco.at(etaBin).at(ptAveBin)->Fill(ptrel,2.*GenJetPt[corrJets(2)]/(GenJetPt[corrJets(0)]+GenJetPt[corrJets(1)]),weight);
	      }
	      if( NobjGenJet > 1 ) {
		double ptGenAve = 0.5*(GenJetColPt[0]+GenJetColPt[1]);
		if( ptGenAve > 0. && NobjGenJet > 2 ) {
		  // Compute imbalance
		  double imbal = 0.;
		  double deltaPhiGen = TVector2::Phi_mpi_pi(GenJetColPhi[0]-GenJetColPhi[1]);
		  double phiDijetAxis = TVector2::Phi_mpi_pi(GenJetColPhi[0]-0.5*deltaPhiGen+M_PI/2.);
		  for(int j = 2; j < NobjGenJet; ++j) {
		    double deltaPhiToDijetAxis = TVector2::Phi_mpi_pi(GenJetColPhi[j]-phiDijetAxis);
		    imbal += cos(deltaPhiToDijetAxis)*GenJetColPt[j];
		  }
		  imbal /= (ptGenAve+std::abs(imbal));
		  
		  //hPt3RelGenVsReco.at(etaBin).at(ptAveBin)->Fill(GenJetColPt[2]/ptGenAve,ptrel,weight);
		  hPt3RelRecoVsImbalGen.at(etaBin).at(ptAveBin)->Fill(imbal,ptrel,weight);
		  hPt3RelGenVsImbalGen.at(etaBin).at(ptAveBin)->Fill(imbal,GenJetColPt[2]/ptGenAve,weight);
		}
	      }
	    }
	  }

	  if( passDeltaPhi && passSoftJet ) {
	    hEta.at(etaBin).at(ptAveBin)->Fill(eta0,weight);
	    hEta.at(etaBin).at(ptAveBin)->Fill(eta1,weight);
	    hPtAve.at(etaBin)->Fill(ptAve,weight);
	    hPtJet.at(etaBin).at(0)->Fill(corrJets.pt(0),weight);
	    hPtJet.at(etaBin).at(1)->Fill(corrJets.pt(1),weight);
	    if( nObjJet > 2 ) {
	      hPtJet.at(etaBin).at(2)->Fill(corrJets.pt(2),weight);
	      if( nObjJet > 3 ) {
		hPtJet.at(etaBin).at(3)->Fill(corrJets.pt(3),weight);
	      }
	    }
	  }

	  if( passDeltaPhi ) {
	    ptAveBin -= binAdmin.hltMinPtBin(hlt,etaBin);
	    newTreesEtaPtAve[etaBin][ptAveBin]->Fill();
	  } else {
	    ++nDeltaPhi;
	  }
	} else {  // End of ptAve bin
	  ++nPtAve;
	}
      } else {	 // End of same eta bin
	++nEta;
      }
    }
     
//     if( !isData ) {
       
//       // GenJet dijet selection
//       if( NobjGenJet > 1 ) {

// 	// Indices of reco-jets matched to gen-jets: the reco-jet closest
// 	// in deltaR is used; this matching has already been done
// 	// when filling the original ntuples
// 	int recoIdx0 = GenJetColJetIdx[0];
// 	int recoIdx1 = GenJetColJetIdx[1];
// 	// Kinematics
// 	double ptGen0 = GenJetColPt[0];
// 	double ptGen1 = GenJetColPt[1];
// 	double ptGenAve = 0.5*(ptGen0+ptGen1);
// 	double ptGen2 = NobjGenJet>2 ? GenJetColPt[2] : 0.;
// 	double pt0 = jetCorrL1[recoIdx0]*jetCorrL2L3[recoIdx0]*jetPt[recoIdx0];
// 	double pt1 = jetCorrL1[recoIdx1]*jetCorrL2L3[recoIdx1]*jetPt[recoIdx1];
// 	// Response of leading two gen-jets
// 	double r0 = 0.;
// 	double r1 = 0.;
// 	if( ptGen0 > 0. && ptGen1 > 0. ) {
// 	  r0 = pt0/ptGen0;
// 	  r1 = pt1/ptGen1;
// 	}
// 	// DeltaR between gen-jet and matched reco-jet
// 	double dr0 = util::deltaR(jetEta[recoIdx0],GenJetColEta[0],jetPhi[recoIdx0],GenJetColPhi[0]);
// 	double dr1 = util::deltaR(jetEta[recoIdx1],GenJetColEta[1],jetPhi[recoIdx1],GenJetColPhi[1]);
// 	bool closeMatch = ( dr0 < maxDeltaR && dr1 < maxDeltaR );
// 	// Do the matched reco-jets pass the jet id?
// 	bool passJetID = ( jetID[recoIdx0] && jetID[recoIdx1] );

// 	if( passJetID && closeMatch ) {
// 	  int ptGenBin = hPtGenBinsForRespVsEtaGen->FindBin(ptGen0)-1;
// 	  if( ptGenBin >= 0 && ptGenBin < static_cast<int>(hRespVsEtaGen.size()) ) {
// 	    double etaGen0 = GenJetColEta[0];
// 	    hRespVsEtaGen.at(ptGenBin)->Fill(etaGen0,r0,wPU);
// 	    hRespVsEtaGenUncorrected.at(ptGenBin)->Fill(etaGen0,jetPt[recoIdx0]/ptGen0,wPU);
// 	    hRespVsEtaGenL1Corrected.at(ptGenBin)->Fill(etaGen0,jetCorrL1[recoIdx0]*jetPt[recoIdx0]/ptGen0,wPU);
// 	  }
// 	  ptGenBin = hPtGenBinsForRespVsEtaGen->FindBin(ptGen1)-1;
// 	  if( ptGenBin >= 0 && ptGenBin < static_cast<int>(hRespVsEtaGen.size()) ) {
// 	    double etaGen1 = GenJetColEta[1];
// 	    hRespVsEtaGen.at(ptGenBin)->Fill(etaGen1,r1,wPU);
// 	    hRespVsEtaGenUncorrected.at(ptGenBin)->Fill(etaGen1,jetPt[recoIdx1]/ptGen1,wPU);
// 	    hRespVsEtaGenL1Corrected.at(ptGenBin)->Fill(etaGen1,jetCorrL1[recoIdx1]*jetPt[recoIdx1]/ptGen1,wPU);
// 	  }
// 	}
	  
// 	unsigned int etaGenBin = 1000;
// 	if( binAdmin.findSameEtaBin(GenJetColEta[0],GenJetColEta[1],etaGenBin) ) {
// 	  double deltaPhi = std::abs(TVector2::Phi_mpi_pi(GenJetColPhi[0]-GenJetColPhi[1]));
// 	  bool passDeltaPhi = deltaPhi > minDeltaPhi;
// 	  int pdgId0 = std::min(21,std::abs(GenPartId_algo[recoIdx0]));
// 	  int pdgId1 = std::min(21,std::abs(GenPartId_algo[recoIdx1]));
	   
// 	  // MC truth response for different selections
// 	  if( passJetID && closeMatch ) {
// 	    // 1) Default
// 	    hRespVsPtGen.at(etaGenBin)->Fill(ptGen0,r0,wPU);
// 	    hRespVsPtGen.at(etaGenBin)->Fill(ptGen1,r1,wPU);
// 	    hRespVsPtGenUncorrected.at(etaGenBin)->Fill(ptGen0,jetPt[recoIdx0]/ptGen0,wPU);
// 	    hRespVsPtGenUncorrected.at(etaGenBin)->Fill(ptGen1,jetPt[recoIdx1]/ptGen1,wPU);
// 	    hRespVsPtGenL1Corrected.at(etaGenBin)->Fill(ptGen0,jetCorrL1[recoIdx0]*jetPt[recoIdx0]/ptGen0,wPU);
// 	    hRespVsPtGenL1Corrected.at(etaGenBin)->Fill(ptGen1,jetCorrL1[recoIdx1]*jetPt[recoIdx1]/ptGen1,wPU);
// 	    hRespVsPtGenVsPdgId.at(etaGenBin).at(pdgId0)->Fill(ptGen0,r0,wPU);
// 	    hRespVsPtGenVsPdgId.at(etaGenBin).at(pdgId1)->Fill(ptGen1,r1,wPU);
// 	    hRespVsPtGenUncorrectedVsPdgId.at(etaGenBin).at(pdgId0)->Fill(ptGen0,jetPt[recoIdx0]/ptGen0,wPU);
// 	    hRespVsPtGenUncorrectedVsPdgId.at(etaGenBin).at(pdgId1)->Fill(ptGen1,jetPt[recoIdx1]/ptGen1,wPU);
// 	    if( PUMCNumVtx < 2 ) {
// 	      hRespVsPtGenUncorrectedLowPU.at(etaGenBin)->Fill(ptGen0,jetPt[recoIdx0]/ptGen0);
// 	      hRespVsPtGenUncorrectedLowPU.at(etaGenBin)->Fill(ptGen1,jetPt[recoIdx1]/ptGen1);
// 	    }

// 	    // 2) Different pile-up scenarios
// 	    hRespVsPtGenPUMC.at(etaGenBin)->Fill(ptGen0,r0);
// 	    hRespVsPtGenPUMC.at(etaGenBin)->Fill(ptGen1,r1);
// 	    if( PUMCNumVtx < 5 ) {
// 	      hRespVsPtGenPULess05.at(etaGenBin)->Fill(ptGen0,r0);
// 	      hRespVsPtGenPULess05.at(etaGenBin)->Fill(ptGen1,r1);
// 	    } else if( PUMCNumVtx < 10 ) {
// 	      hRespVsPtGenPULess10.at(etaGenBin)->Fill(ptGen0,r0);
// 	      hRespVsPtGenPULess10.at(etaGenBin)->Fill(ptGen1,r1);
// 	    } else if( PUMCNumVtx < 15 ) {
// 	      hRespVsPtGenPULess15.at(etaGenBin)->Fill(ptGen0,r0);
// 	      hRespVsPtGenPULess15.at(etaGenBin)->Fill(ptGen1,r1);
// 	    } else {
// 	      hRespVsPtGenPULess99.at(etaGenBin)->Fill(ptGen0,r0);
// 	      hRespVsPtGenPULess99.at(etaGenBin)->Fill(ptGen1,r1);
// 	    }
// 	  }

// 	  // 3) JetID
// 	  if( closeMatch ) {
// 	    hRespVsPtGenNoJetID.at(etaGenBin)->Fill(ptGen0,r0,wPU);
// 	    hRespVsPtGenNoJetID.at(etaGenBin)->Fill(ptGen1,r1,wPU);
// 	  }

// 	  // 4) DeltaR matching
// 	  if( passJetID ) {
// 	    hDeltaRVsPtGen.at(etaGenBin)->Fill(ptGen0,dr0);
// 	    hDeltaRVsPtGen.at(etaGenBin)->Fill(ptGen1,dr1);
	     
// 	    hRespVsPtGenNoDeltaR.at(etaGenBin)->Fill(ptGen0,r0,wPU);
// 	    hRespVsPtGenNoDeltaR.at(etaGenBin)->Fill(ptGen1,r1,wPU);

// 	    if( dr0 < 0.3 && dr1 < 0.3 ) {
// 	      hRespVsPtGenDeltaRLess30.at(etaGenBin)->Fill(ptGen0,r0,wPU);
// 	      hRespVsPtGenDeltaRLess30.at(etaGenBin)->Fill(ptGen1,r1,wPU);
// 	      if( dr0 < 0.25 && dr1 < 0.25 ) {
// 		hRespVsPtGenDeltaRLess25.at(etaGenBin)->Fill(ptGen0,r0,wPU);
// 		hRespVsPtGenDeltaRLess25.at(etaGenBin)->Fill(ptGen1,r1,wPU);
// 		if( dr0 < 0.2 && dr1 < 0.2 ) {
// 		  hRespVsPtGenDeltaRLess20.at(etaGenBin)->Fill(ptGen0,r0,wPU);
// 		  hRespVsPtGenDeltaRLess20.at(etaGenBin)->Fill(ptGen1,r1,wPU);
// 		  if( dr0 < 0.15 && dr1 < 0.15 ) {
// 		    hRespVsPtGenDeltaRLess15.at(etaGenBin)->Fill(ptGen0,r0,wPU);
// 		    hRespVsPtGenDeltaRLess15.at(etaGenBin)->Fill(ptGen1,r1,wPU);
// 		    if( dr0 < 0.1 && dr1 < 0.1 ) {
// 		      hRespVsPtGenDeltaRLess10.at(etaGenBin)->Fill(ptGen0,r0,wPU);
// 		      hRespVsPtGenDeltaRLess10.at(etaGenBin)->Fill(ptGen1,r1,wPU);
// 		      if( dr0 < 0.05 && dr1 < 0.05 ) {
// 			hRespVsPtGenDeltaRLess05.at(etaGenBin)->Fill(ptGen0,r0,wPU);
// 			hRespVsPtGenDeltaRLess05.at(etaGenBin)->Fill(ptGen1,r1,wPU);
// 		      }
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }

// 	  // 5) Dijet selection
// 	  if( passDeltaPhi ) { 
// 	    // Pt soft cuts
// 	    for(int ptSoftBin = binAdmin.nPtSoftBins()-1; ptSoftBin >= 0; --ptSoftBin) {
// 	      if( ptGen2 > binAdmin.ptSoftMax(ptSoftBin)*ptGenAve ) break;
// 	      if( passJetID && closeMatch ) {
// 		hRespVsPtGenVsPtSoft.at(etaGenBin).at(ptSoftBin)->Fill(ptGen0,r0,wPU);
// 		hRespVsPtGenVsPtSoft.at(etaGenBin).at(ptSoftBin)->Fill(ptGen1,r1,wPU);
// 	      }
// 	      // Gen-jet spectrum independent from reco-jet properties
// 	      hPtGen.at(etaGenBin).at(ptSoftBin)->Fill(ptGen0,wSpec);
// 	      hPtGen.at(etaGenBin).at(ptSoftBin)->Fill(ptGen1,wSpec);
// 	      hPtGenUnweighted.at(etaGenBin).at(ptSoftBin)->Fill(ptGen0);
// 	      hPtGenUnweighted.at(etaGenBin).at(ptSoftBin)->Fill(ptGen1);
// 	    }
// 	    if( passJetID && closeMatch && ptGen2 < 0.14*ptGenAve ) {
// 	      hRespVsPtGenVsPdgIdDijets.at(etaGenBin).at(pdgId0)->Fill(ptGen0,r0,wPU);
// 	      hRespVsPtGenVsPdgIdDijets.at(etaGenBin).at(pdgId1)->Fill(ptGen1,r1,wPU);
// 	    }
	      
// 	    if( writeAllTrees ) {
// 	      // Find ptGenAve bin
// 	      unsigned int ptGenAveBin = 1000;
// 	      if( binAdmin.findPtBin(ptGenAve,etaGenBin,ptGenAveBin) ) {
// 		for(int j = 0; j < nObjJet; ++j) {
// 		  corrJetIdx[j] = GenJetColJetIdx[j];
// 		}
// 		newTreesEtaGenPtGenAve[etaGenBin][ptGenAveBin]->Fill();
// 	      }
// 	    }
// 	  } // End of if( sameEtaBin )
// 	}
//   }
// }
   
     
    if( i%50000 == 0 ) {
      for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
  	for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(hlt,etaBin); ++ptBin) {
 	  newTreesEtaPtAve[etaBin][ptBin]->AutoSave();
 	  if( writeAllTrees ) {
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
	newTreesEtaGenPtGenAve[etaBin][ptBin]->AutoSave();
      }
    }
  }

  TString fileControlPlotsName = outFilePath+"/Kalibri_ControlPlots";
  if( maxHltThres > 0 ) {
    fileControlPlotsName += "_HLTDiJetAve"+util::toTString(maxHltThres);
  }
  fileControlPlotsName += ".root";
  TFile fileControlPlots(fileControlPlotsName,"RECREATE");
  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
      fileControlPlots.WriteTObject(hEta.at(etaBin).at(ptBin));
      fileControlPlots.WriteTObject(hEtaNMin1.at(etaBin).at(ptBin));
      fileControlPlots.WriteTObject(hDeltaPhiNMin1.at(etaBin).at(ptBin));
      fileControlPlots.WriteTObject(hPtRelNMin1.at(etaBin).at(ptBin));
      fileControlPlots.WriteTObject(hDeltaPhiVsPt3Rel.at(etaBin).at(ptBin));
      if( !isData ) {
	fileControlPlots.WriteTObject(hPt3RelGenVsReco.at(etaBin).at(ptBin));
	fileControlPlots.WriteTObject(hPt3RelRecoVsImbalGen.at(etaBin).at(ptBin));
	fileControlPlots.WriteTObject(hPt3RelGenVsImbalGen.at(etaBin).at(ptBin));
      }
      fileControlPlots.WriteTObject(hPtAve.at(etaBin));
      for(unsigned int i = 0; i < hPtJet.at(etaBin).size(); ++i) {
	fileControlPlots.WriteTObject(hPtJet.at(etaBin).at(i));
      }
    }
  }

  if( !isData ) {
    TFile fileResponse(outFilePath+"/Kalibri_MCTruthResponse.root","RECREATE");
    TFile fileSpetrum(outFilePath+"/Kalibri_MCTruthSpectrum.root","RECREATE");
    fileResponse.WriteTObject(hPtGenBinsForRespVsEtaGen);
    for(unsigned int ptBin = 0; ptBin < hRespVsEtaGen.size(); ++ptBin) {
      fileResponse.WriteTObject(hRespVsEtaGen.at(ptBin));
      fileResponse.WriteTObject(hRespVsEtaGenUncorrected.at(ptBin));
      fileResponse.WriteTObject(hRespVsEtaGenL1Corrected.at(ptBin));
    }
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      fileResponse.WriteTObject(hRespVsPtGen.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenUncorrected.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenUncorrectedLowPU.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenL1Corrected.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenPUMC.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenPULess05.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenPULess10.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenPULess15.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenPULess99.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenNoJetID.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenNoDeltaR.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenDeltaRLess05.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenDeltaRLess10.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenDeltaRLess15.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenDeltaRLess20.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenDeltaRLess25.at(etaBin));
      fileResponse.WriteTObject(hRespVsPtGenDeltaRLess30.at(etaBin));
      fileResponse.WriteTObject(hDeltaRVsPtGen.at(etaBin));
      
      for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
	fileResponse.WriteTObject(hRespVsPtGenVsPtSoft.at(etaBin).at(ptSoftBin));
	fileSpetrum.WriteTObject(hPtGen.at(etaBin).at(ptSoftBin));
	fileSpetrum.WriteTObject(hPtGenUnweighted.at(etaBin).at(ptSoftBin));
      }
      for(unsigned int id = 0; id < hRespVsPtGenVsPdgId.at(etaBin).size(); ++id) {
	fileResponse.WriteTObject(hRespVsPtGenVsPdgId.at(etaBin).at(id));
	fileResponse.WriteTObject(hRespVsPtGenUncorrectedVsPdgId.at(etaBin).at(id));
	fileResponse.WriteTObject(hRespVsPtGenVsPdgIdDijets.at(etaBin).at(id));
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
  std::cout << "  " << (nEntries -= nJetID ) << " events with Jet(1,2) passing loose JetID cuts" << std::endl;
  if( !isData ) std::cout << "  " << (nEntries -= nLargeResp ) << " events with good response" << std::endl;
  std::cout << "  " << (nEntries -= nDeltaPhi ) << " events with |DeltaPhi(1,2)| > " << minDeltaPhi << std::endl;
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
	newFilesEtaGenPtGenAve[etaBin][ptBin]->Close();
	delete newFilesEtaGenPtGenAve[etaBin][ptBin];
      }
    }
  }
  delete oldChain;

  for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
    for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
      delete hEta.at(etaBin).at(ptBin);
      delete hEtaNMin1.at(etaBin).at(ptBin);
      delete hDeltaPhiNMin1.at(etaBin).at(ptBin);
      delete hPtRelNMin1.at(etaBin).at(ptBin);
      delete hDeltaPhiVsPt3Rel.at(etaBin).at(ptBin);
      delete hPt3RelGenVsReco.at(etaBin).at(ptBin);
      delete hPt3RelRecoVsImbalGen.at(etaBin).at(ptBin);
      delete hPt3RelGenVsImbalGen.at(etaBin).at(ptBin);
    }
    delete hPtAve.at(etaBin);
    for(unsigned int i = 0; i < hPtJet.at(etaBin).size(); ++i) {
      delete hPtJet.at(etaBin).at(i);
    }
  }

  if( !isData ) {
    delete hPtGenBinsForRespVsEtaGen;
    for(unsigned int ptBin = 0; ptBin < hRespVsEtaGen.size(); ++ptBin) {
      delete hRespVsEtaGen.at(ptBin);
      delete hRespVsEtaGenUncorrected.at(ptBin);
      delete hRespVsEtaGenL1Corrected.at(ptBin);
    }
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      delete hRespVsPtGen.at(etaBin);
      delete hRespVsPtGenUncorrected.at(etaBin);
      delete hRespVsPtGenUncorrectedLowPU.at(etaBin);
      delete hRespVsPtGenL1Corrected.at(etaBin);
      delete hRespVsPtGenPUMC.at(etaBin);
      delete hRespVsPtGenPULess05.at(etaBin);
      delete hRespVsPtGenPULess10.at(etaBin);
      delete hRespVsPtGenPULess15.at(etaBin);
      delete hRespVsPtGenPULess99.at(etaBin);
      delete hRespVsPtGenNoJetID.at(etaBin);
      delete hRespVsPtGenNoDeltaR.at(etaBin);
      delete hRespVsPtGenDeltaRLess05.at(etaBin);
      delete hRespVsPtGenDeltaRLess10.at(etaBin);
      delete hRespVsPtGenDeltaRLess15.at(etaBin);
      delete hRespVsPtGenDeltaRLess20.at(etaBin);
      delete hRespVsPtGenDeltaRLess25.at(etaBin);
      delete hRespVsPtGenDeltaRLess30.at(etaBin);
      delete hDeltaRVsPtGen.at(etaBin);
  
      for(unsigned int ptSoftBin = 0; ptSoftBin < binAdmin.nPtSoftBins(); ++ptSoftBin) {
	delete hRespVsPtGenVsPtSoft.at(etaBin).at(ptSoftBin);
	delete hPtGen.at(etaBin).at(ptSoftBin);
	delete hPtGenUnweighted.at(etaBin).at(ptSoftBin);
      }
      for(unsigned int id = 0; id < hRespVsPtGenVsPdgId.at(etaBin).size(); ++id) {
	delete hRespVsPtGenVsPdgId.at(etaBin).at(id);
	delete hRespVsPtGenUncorrectedVsPdgId.at(etaBin).at(id);
	delete hRespVsPtGenVsPdgIdDijets.at(etaBin).at(id);
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
