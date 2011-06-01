// $Id:  $


#include <iostream>
#include <vector>

#include "TBranch.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include "BinningAdmin.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/FileOps.h"



void reweightEventsToPUDistribution(const char* path, const TString &fileNameDataPUDist);
void reweightEventsToMCStats(const char* path);
std::vector<double> generate_flat10_weights(const TH1* data_npu_estimated);



// Reweight events to reflect the data PU distribution
// --------------------------------------------------
void reweightEventsToPUDistribution(const char* path, const TString &fileNameDataPUDist) {


  // ++++ Set up +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Setup" << std::endl;

  std::cout << "  Computing PU weights" << std::endl;
  TH1* hDataPU = util::FileOps::readTH1(fileNameDataPUDist,"pileup");
  std::vector<double> weightsPU = generate_flat10_weights(hDataPU);

  std::cout << "  Adding files to TChain" << std::endl;
  TChain* c = new TChain("DiJetTree","DiJetTree");
  c->Add(path);
      


  // ++++ Reweight +++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TObjArray* fileElements = c->GetListOfFiles();
  TIter next(fileElements);
  TChainElement* chEl = 0;
  while( ( chEl=(TChainElement*)next() ) ) {
    std::cout << "Processing file '" << chEl->GetTitle() << "'\n";

    std::cout << "  Preparing tree" << std::endl;
    // Get tree from file
    TFile f(chEl->GetTitle(),"update");
    TTree* djt = (TTree*)f.Get("DiJetTree");

    Float_t oldWeight = 0.;
    Float_t newWeight = 0.;
    Int_t PUMCNumVtx = 0;
    djt->SetBranchAddress("Weight",&oldWeight);
    djt->SetBranchAddress("PUMCNumVtx",&PUMCNumVtx);

    // Clone tree and set new weight branch
    djt->SetBranchStatus("Weight",0);
    TTree* ndjt = djt->CloneTree(-1,"fast");
    TBranch* bWeight = ndjt->Branch("Weight",&newWeight,"Weight/F");
    djt->SetBranchStatus("Weight",1);

    // Set new weight
    std::cout << "  Setting new weights" << std::endl;
    for(Long64_t i = 0; i < djt->GetEntries(); ++i) {
      djt->GetEntry(i);
      if( PUMCNumVtx < static_cast<int>(weightsPU.size()) ) {
	newWeight = oldWeight * weightsPU.at(PUMCNumVtx);
      } else {
	std::cerr << "ERROR: no event weights for PUMCNumVtx = " << PUMCNumVtx << std::endl;
	newWeight = oldWeight;
      }
      bWeight->Fill();
    }

    std::cout << "  Cleaning up" << std::endl;
    // Delete original tree and replace by new one
    delete djt;
    f.Delete("DiJetTree;1");
    ndjt->Write("", TObject::kOverwrite);
    f.Close();
  }
  delete c;
}




// Reweight events such that the largest weight is one
// --------------------------------------------------
void reweightEventsToMCStats(const char* path) {


  // ++++ Set up +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  std::cout << "Setup" << std::endl;
  TChain* c = new TChain("DiJetTree","DiJetTree");
  c->Add(path);
  const double maxResp = 2.5;
      


  // ++++ Reweight +++++++++++++++++++++++++++++++++++++++++++++++++++++++

  TObjArray* fileElements = c->GetListOfFiles();
  TIter next(fileElements);
  TChainElement* chEl = 0;
  while( ( chEl=(TChainElement*)next() ) ) {
    std::cout << "Processing file '" << chEl->GetTitle() << "'\n";

    // Get tree from file
    TFile f(chEl->GetTitle(),"update");
    TTree* djt = (TTree*)f.Get("DiJetTree");

    // Find largest weight
    std::cout << "  Searching for largest weight (in events with r < " << maxResp << ")" << std::endl;
    Float_t largestWeight = 0.;
    Float_t oldWeight = 0.;
    Float_t newWeight = 0.;
    Int_t NObjJet = 0;
    const int maxNJet = 50;
    Int_t Idx[maxNJet];
    Float_t JetPt[maxNJet];
    Float_t GenJetPt[maxNJet];
    Float_t JetCorrL1[maxNJet];
    Float_t JetCorrL2L3[maxNJet];
    djt->SetBranchAddress("NobjJet",&NObjJet);
    djt->SetBranchAddress("JetPt",JetPt);
    djt->SetBranchAddress("GenJetPt", GenJetPt);
    djt->SetBranchAddress("JetCorrL1",JetCorrL1);
    djt->SetBranchAddress("JetCorrL2L3",JetCorrL2L3);
    djt->SetBranchAddress("L2L3CorrJetColJetIdx",Idx);
    djt->SetBranchAddress("Weight", &oldWeight);
    Long64_t nentries = djt->GetEntries();
    for(Long64_t i = 0; i < nentries; i++){
      djt->GetEntry(i);
      double p1 = JetCorrL1[Idx[0]]*JetCorrL2L3[Idx[0]]*JetPt[Idx[0]];
      double p2 = JetCorrL1[Idx[1]]*JetCorrL2L3[Idx[1]]*JetPt[Idx[1]];
      double p3 = NObjJet > 2 ? JetCorrL1[Idx[2]]*JetCorrL2L3[Idx[2]]*JetPt[Idx[2]] : 0.;
      double pg1 = GenJetPt[Idx[0]];
      double pg2 = GenJetPt[Idx[1]];
      double r1 = 0.;
      double r2 = 0.;
      if( pg1 && pg2 && p3 < 0.225*0.5*(p1+p2) ) {
	r1 = p1/pg1;
	r2 = p2/pg2;
	if( (r1 < maxResp) && (r2 < maxResp) && (oldWeight > largestWeight) ) largestWeight = oldWeight;
      }
    }
    if( !(largestWeight > 0.) ) largestWeight = 1.;
    std::cout << "    Largest weight = " << largestWeight << std::endl;
    
    // Clone tree and set new weight branch
    djt->SetBranchStatus("Weight",0);
    TTree* ndjt = djt->CloneTree(-1,"fast");
    TBranch* bWeight = ndjt->Branch("Weight", &newWeight, "Weight/F");
    djt->SetBranchStatus("Weight",1);

    // Set new weight
    std::cout << "  Setting new weights 'w --> w/" << largestWeight << "'" << std::endl;
    for (Long64_t i = 0; i < nentries; i++){
      djt->GetEntry(i);
      newWeight = oldWeight/largestWeight;
      bWeight->Fill();
    }

    // Delete original tree and replace by new one
    delete djt;
    f.Delete("DiJetTree;1");
    ndjt->Write("", TObject::kOverwrite);
    f.Close();
  }
  delete c;
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
