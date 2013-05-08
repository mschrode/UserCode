// $Id: $
//

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
#include "TRandom.h"
#include "TRandom1.h"
#include "TROOT.h"
#include "TString.h"
#include "TVector2.h"

#include "BinningAdmin.h"
#include "WeightProducer.h"
#include "Parameters.h"
#define UTILS_AS_HEADER_FILE
#include "../util/utils.h"
#include "../util/FileOps.h"
#include "../util/HistOps.h"



namespace sampleTools {

  // ++++ PARAMETER DEFINITIONS ++++++++++++++++++++++++++++++++

  // Event selection cuts
  const double minDeltaPhi_ = 2.7;
  const double maxDeltaR_ = 0.1;




  // ++++ FUNCTION DELCARATIONS ++++++++++++++++++++++++++++++++

  // ++++ Functions ++++++++++++++++++++++++++++++++++++++++++++
  TChain *createTChain(const TString &fileListName);
  void processEvents();



  // ++++ MAIN FUNCTION ++++++++++++++++++++++++++++++++++++++++
  void fillResTailInputHists() {
    if( Parameters::isInit ) {
      processEvents();
    } else {
      std::cerr << "ERROR: Parameters not initialized" << std::endl;
      exit(-1);
    }
  }




  // ++++ FUNCTION IMPLEMENTATIONS +++++++++++++++++++++++++++++

  // --------------------------------------------------
  void processEvents() {
  
    // ++++ Set parameters +++++++++++++++++++++++++++++++++++++++
  
    std::cout << "Initializing parameters" << std::endl;

    const TString inFileList = Parameters::runOnData ? Parameters::inFileListData : Parameters::inFileListMC;

    WeightProducer::init(Parameters::qcdSampleType,Parameters::puScenario,util::FileOps::readTH1(Parameters::puHistFile,"pileup"));

    // Cut-flow counters
    unsigned int nMaxNJet = 0;
    unsigned int nDijets = 0;
    unsigned int nHlt = 0;
    unsigned int nJetID = 0;
    unsigned int nLargeResp = 0;
    unsigned int nPt3Rel = 0;
    unsigned int nDeltaPhi = 0;
    unsigned int nEta = 0;
    unsigned int nPtAve = 0;
  
    // Event weight
    double weight = 1.;   // Combined weight
    double wPU = 1.;	// PU weight
    float wSpec = 1.;	// Spectrum weight (for flat samples, taken from input ntuple)

    TRandom* rand_ = new TRandom1(0);
  
  

    // ++++ Prepare objects ++++++++++++++++++++++++++++++++++++++

    std::cout << "Defining binning and trigger thresholds" << std::endl;

    const sampleTools::BinningAdmin binAdmin(Parameters::configAdm,Parameters::configBin);
    binAdmin.printBinning();
    if( Parameters::runOnData ) binAdmin.print();


    std::cout << "Preparing histograms" << std::endl;

    // Asymmetry distributions
    std::vector< std::vector< std::vector< TH1* > > > hPtAbsAsym;	// [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector< std::vector< TH1* > > > hPtAsym; // [pt3RelBin][etaBin][ptAveBin]

    // Control plots
    std::vector< std::vector< std::vector<TH1*> > > hPtAve; // [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector< std::vector<TH1*> > > hResp; // [pt3RelBin][etaBin][ptAveBin]


    std::vector< std::vector<TH1*> > hPtJet; // [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector<TH1*> > hEta;
    std::vector< std::vector<TH1*> > hEtaNMin1;
    std::vector< std::vector<TH1*> > hDeltaPhiNMin1;
    std::vector< std::vector<TH1*> > hPt3RelNMin1;
    std::vector< std::vector<TH1*> > hPtGen; // [pt3RelBin][etaGenBin]
  
    // MC-truth response and control plots
    std::vector< std::vector<TH2*> > hRespVsPtGen; // [pt3RelBin][etaGenBin]
    std::vector< std::vector<TH2*> > hRespVsEtaGen; // [pt3RelBin][ptGenBin]
    std::vector< std::vector<TH2*> > hRespVsPtGenPUMC;
    std::vector< std::vector<TH2*> > hRespVsPtGenPULess05;
    std::vector< std::vector<TH2*> > hRespVsPtGenPULess10;
    std::vector< std::vector<TH2*> > hRespVsPtGenPULess15;
    std::vector< std::vector<TH2*> > hRespVsPtGenPULess99;
    std::vector< std::vector<TH2*> > hRespVsPtGenDeltaRLess05;
    std::vector< std::vector<TH2*> > hRespVsPtGenDeltaRLess10;
    std::vector< std::vector<TH2*> > hRespVsPtGenDeltaRLess15;
    std::vector< std::vector<TH2*> > hRespVsPtGenDeltaRLess20;
    std::vector< std::vector<TH2*> > hRespVsPtGenDeltaRLess25;
    std::vector< std::vector<TH2*> > hRespVsPtGenDeltaRLess30;
    std::vector< std::vector<TH2*> > hDeltaRVsPtGen;


    // Equidistant log bins
    std::vector<double> ptBinEdges(26);
    util::HistOps::equidistLogBins(ptBinEdges,ptBinEdges.size()-1,4.,2700.);

    // Loop over pt3RelBins
    for(unsigned int pt3RelBin = 0; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {

      std::vector< std::vector< TH1* > > hPtAbsAsymTmp2D;
      std::vector< std::vector< TH1* > > hPtAsymTmp2D;
      std::vector< std::vector< TH1* > > hPtAveTmp2D;
      std::vector< std::vector< TH1* > > hRespTmp2D;

      // Loop over etaBins
      for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {

	std::vector< TH1* > hPtAbsAsymTmp;
	std::vector< TH1* > hPtAsymTmp;
	std::vector< TH1* > hPtAveTmp;
	std::vector< TH1* > hRespTmp;

	// Loop over ptAveBins
	for(unsigned int ptAveBin = 0; ptAveBin < binAdmin.nPtBins(etaBin); ++ptAveBin) {

	  // Plots per (pt3Rel,eta,ptAve) bin
	  TString binId = 
	    "_Eta"+util::toTString(etaBin) +
	    "_Pt"+util::toTString(ptAveBin) +
	    "_PtSoft"+util::toTString(pt3RelBin);

	  TString title = 
	    util::toTString(binAdmin.etaMin(etaBin)) + 
	    " < |#eta| < " + 
	    util::toTString(binAdmin.etaMax(etaBin)) +
	    ",  " +
	    util::toTString(binAdmin.ptMin(etaBin,ptAveBin)) +
	    " < p^{ave}_{T} < " +
	    util::toTString(binAdmin.ptMax(etaBin,ptAveBin)) +	
	    " GeV,  " + util::LabelFactory::pt3RelMax() + " < " +
	    util::toTString(binAdmin.ptSoftMax(pt3RelBin));

	  TH1* h = new TH1D(Parameters::hNamePtAsymAbs+binId,title+";|p_{T} asymmetry|;Events",150,-1.,1.);
	  h->Sumw2();
	  hPtAbsAsymTmp.push_back(h);

	  h = new TH1D(Parameters::hNamePtAsym+binId,title+";p_{T} asymmetry;Events",101,-1.,1.);
	  h->Sumw2();
	  hPtAsymTmp.push_back(h);

	  h = new TH1D(Parameters::hNamePtAve+binId,title+";p^{ave}_{T} (GeV);Events",800,0,2000);
	  h->Sumw2();
	  hPtAveTmp.push_back(h);

	  h = new TH1D(Parameters::hNameResp+binId,title+";Response;Jets",120,0,2);
	  h->Sumw2();
	  hRespTmp.push_back(h);

	}	// End of loop over ptAveBins

	hPtAbsAsymTmp2D.push_back(hPtAbsAsymTmp);
	hPtAsymTmp2D.push_back(hPtAsymTmp);
	hPtAveTmp2D.push_back(hPtAveTmp);
	hRespTmp2D.push_back(hRespTmp);

      } // End of loop over etaBin

      hPtAbsAsym.push_back(hPtAbsAsymTmp2D);
      hPtAsym.push_back(hPtAsymTmp2D);
      hPtAve.push_back(hPtAveTmp2D);
      hResp.push_back(hRespTmp2D);

    } // End of loop over pt3RelBins



    std::cout << "Preparing trees" << std::endl;

    // Open Kalibri ntuples
    TChain *chain = createTChain(inFileList);
  
    // Deactivate branches not needed
    chain->SetBranchStatus("Track*",0);
    chain->SetBranchStatus("NobjTrack",0);
    chain->SetBranchStatus("Tow*",0);
    chain->SetBranchStatus("GenPart*",0);
    chain->SetBranchStatus("GenPartId*",1);
    chain->SetBranchStatus("Vtx*",0);
    chain->SetBranchStatus("VtxN",1);
    chain->SetBranchStatus("Met*",0);
    chain->SetBranchStatus("Mu*",0);

    // Declare elements
    const int maxNJet = 50;
    Int_t NObjJet = 0;
    Float_t JetPt[maxNJet];
    Float_t JetEta[maxNJet];
    Float_t JetPhi[maxNJet];
    Float_t JetCorrL1[maxNJet];
    Float_t JetCorrL2L3[maxNJet];	// Should contain L2*L3*LResidual, starting from L1 corrected pt
    Float_t JetCorrUncert[maxNJet];
    Bool_t JetID[maxNJet];
    Bool_t HltDiPFJetAve40 = false;
    Bool_t HltDiPFJetAve80 = false;
    Bool_t HltDiPFJetAve140 = false;
    Bool_t HltDiPFJetAve200 = false;
    Bool_t HltDiPFJetAve260 = false;
    Bool_t HltDiPFJetAve320 = false;
    Bool_t HltDiPFJetAve400 = false;

    //   // Get other elements
    //   UInt_t          RunNumber;
    //   UInt_t          LumiBlockNumber;
    //   UInt_t          EventNumber;
    //   Int_t           NobjTow;
    //   Float_t         JetEt[maxNJet];   //[NobjJet]
    //   Float_t         JetE[maxNJet];   //[NobjJet]
    //   Int_t           JetN90Hits[maxNJet];   //[NobjJet]
    //   Float_t         JetHad[maxNJet];   //[NobjJet]
    //   Float_t         JetEMF[maxNJet];   //[NobjJet]
    //   Float_t         JetFHPD[maxNJet];   //[NobjJet]
    //   Float_t         JetFRBX[maxNJet];   //[NobjJet]
    //   Float_t         JetEtWeightedSigmaPhi[maxNJet];   //[NobjJet]
    //   Float_t         JetEtWeightedSigmaEta[maxNJet];   //[NobjJet]
    //   Float_t         JetCorrZSP[maxNJet];   //[NobjJet]
    //   Float_t         JetCorrL2[maxNJet];   //[NobjJet]
    //   Float_t         JetCorrL3[maxNJet];   //[NobjJet]
    //   Float_t         JetCorrJPT[maxNJet];   //[NobjJet]
    //   Float_t         JetCorrL2L3JPT[maxNJet];   //[NobjJet]
    //   Float_t         JetCorrL4JW[maxNJet];   //[NobjJet]
    //   Int_t           JetIEta[maxNJet];   //[NobjJet]
    //   Int_t           JetIPhi[maxNJet];   //[NobjJet]
    //   Float_t         JetGenJetDeltaR[maxNJet];   //[NobjJet]
    //   Float_t         GenJetEt[maxNJet];   //[NobjJet]
    //   Float_t         GenJetE[maxNJet];   //[NobjJet]

    Float_t         GenJetPt[maxNJet];   //[NobjJet]
    Float_t         GenJetPhi[maxNJet];   //[NobjJet]
    Float_t         GenJetEta[maxNJet];   //[NobjJet]

    Int_t           NobjGenJet;
    Float_t         GenJetColPt[maxNJet];   //[NobjGenJet]
    Float_t         GenJetColPhi[maxNJet];   //[NobjGenJet]
    Float_t         GenJetColEta[maxNJet];   //[NobjGenJet]
    //   Float_t         GenJetColEt[maxNJet];   //[NobjGenJet]
    //   Float_t         GenJetColE[maxNJet];   //[NobjGenJet]
    //   Float_t         GenJetColEmE[maxNJet];   //[NobjGenJet]
    //   Float_t         GenJetColHadE[maxNJet];   //[NobjGenJet]
    //   Float_t         GenJetColInvE[maxNJet];   //[NobjGenJet]
    //   Float_t         GenJetColAuxE[maxNJet];   //[NobjGenJet]
    Int_t           GenJetColJetIdx[maxNJet];   //[NobjGenJet]
    Float_t         GenEvtScale;
    //   Int_t           VtxN = 0;
    Int_t           PUMCNumVtx = 0;
    //   Int_t           GenPartId_algo[maxNJet];   //[NobjJet]
    //   Int_t           GenPartId_phys[maxNJet];   //[NobjJet]


    // Set branch addresses
    chain->SetBranchAddress("NobjJet",&NObjJet);
    chain->SetBranchAddress("JetPt",JetPt);
    chain->SetBranchAddress("JetEta",JetEta);
    chain->SetBranchAddress("JetPhi",JetPhi);
    chain->SetBranchAddress("JetCorrL1",JetCorrL1);
    chain->SetBranchAddress("JetCorrL2L3",JetCorrL2L3);
    chain->SetBranchAddress("HltDiPFJetAve40",&HltDiPFJetAve40);
    chain->SetBranchAddress("HltDiPFJetAve80",&HltDiPFJetAve80);
    chain->SetBranchAddress("HltDiPFJetAve140",&HltDiPFJetAve140);
    chain->SetBranchAddress("HltDiPFJetAve200",&HltDiPFJetAve200);
    chain->SetBranchAddress("HltDiPFJetAve260",&HltDiPFJetAve260);
    chain->SetBranchAddress("HltDiPFJetAve320",&HltDiPFJetAve320);
    chain->SetBranchAddress("HltDiPFJetAve400",&HltDiPFJetAve400);
    chain->SetBranchAddress("JetIDLoose",JetID);
    chain->SetBranchAddress("JetCorrUncert",JetCorrUncert);
    chain->SetBranchAddress("NobjGenJet", &NobjGenJet);
    chain->SetBranchAddress("GenJetPt", GenJetPt);
    chain->SetBranchAddress("GenJetPhi", GenJetPhi);
    chain->SetBranchAddress("GenJetEta", GenJetEta);
    chain->SetBranchAddress("GenJetColPt", GenJetColPt);
    chain->SetBranchAddress("GenJetColPhi", GenJetColPhi);
    chain->SetBranchAddress("GenJetColEta", GenJetColEta);
    chain->SetBranchAddress("GenJetColJetIdx", GenJetColJetIdx);
    chain->SetBranchAddress("GenEvtScale",&GenEvtScale);
    chain->SetBranchAddress("PUMCNumVtx",&PUMCNumVtx);



    // Container for jet ordering
    util::JetIndexCol corrJets;


    // Container for trigger decisions
    // For each eta-pt bin, consider only decision of highest,
    // unprescaled trigger
    std::vector< std::vector<bool*> > hltDecisions(binAdmin.nEtaBins());
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      std::vector<bool*> hltDecisionsInEtaBin(binAdmin.nPtBins(etaBin));
      for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
	TString hltName = "";
	if( binAdmin.findHltMax(etaBin,ptBin,hltName) ) {	// Knows the trigger turn-on threshold
	
	  if( hltName == "HLT_DiPFJetAve40" ) 
	    hltDecisionsInEtaBin.at(ptBin) = &HltDiPFJetAve40;   
	  else if( hltName == "HLT_DiPFJetAve80" ) 
	    hltDecisionsInEtaBin.at(ptBin) = &HltDiPFJetAve80;   
	  else if( hltName == "HLT_DiPFJetAve140" ) 
	    hltDecisionsInEtaBin.at(ptBin) = &HltDiPFJetAve140;   
	  else if( hltName == "HLT_DiPFJetAve200" ) 
	    hltDecisionsInEtaBin.at(ptBin) = &HltDiPFJetAve200;   
	  else if( hltName == "HLT_DiPFJetAve260" ) 
	    hltDecisionsInEtaBin.at(ptBin) = &HltDiPFJetAve260;   
	  else if( hltName == "HLT_DiPFJetAve320" ) 
	    hltDecisionsInEtaBin.at(ptBin) = &HltDiPFJetAve320;   
	  else if( hltName == "HLT_DiPFJetAve400" ) 
	    hltDecisionsInEtaBin.at(ptBin) = &HltDiPFJetAve400;
	  else {
	    std::cerr << "ERROR setting up trigger-bin matching" << std::endl;
	    std::cerr << "  HLT '" << hltName << "' not found" << std::endl;
	    exit(-1);
	  }
	
	} else {
	  std::cerr << "ERROR setting up trigger-bin matching" << std::endl;
	  std::cerr << "  HLT turn-on does not match eta-pt bin" << std::endl;
	  exit(-1);
	}
      }
      hltDecisions.at(etaBin) = hltDecisionsInEtaBin;
    }
  



    // ++++ Loop over tree and select dijets +++++++++++++++++

    std::cout << "Processing events" << std::endl;
    int nEntries = chain->GetEntries();
    if( Parameters::nEvts > 0 && Parameters::nEvts <= nEntries ) nEntries = Parameters::nEvts;

    for(int i = 0; i < nEntries; ++i) {
      if( i%50000 == 0 ) {
	std::cout << "  Processed " << i << " events" << std::endl;
      }

      chain->GetEntry(i);

      if( NObjJet > maxNJet ) {
	std::cerr << "WARNING: nObjJet = " << NObjJet << " > " << maxNJet << ". Skipping event!\n";
	++nMaxNJet;
	continue;
      }


      // Weights
      if( !Parameters::runOnData ) {
	// Spectrum re-weighting
	wSpec = WeightProducer::qcdSpectrumWeight(GenEvtScale);
	// PU re-weighting
	wPU = WeightProducer::puWeight(PUMCNumVtx);
	// Combined event weight: spectrum weight x PU weight
	weight = wSpec * wPU;
      }


      // Order corrected jets
      corrJets.clear();
      for(int j = 0; j < NObjJet; ++j) {
	corrJets.add(j,JetCorrL1[j]*JetCorrL2L3[j]*JetPt[j]);
      }
      corrJets.sort();
     

      // Dijet selection
      if( NObjJet < 2 ) {
	++nDijets;
	continue;
      }


      // Reject MC events with large response
      // These events were introduced in 2011 MC due to incorrect
      // PU mixing
      if ( !Parameters::runOnData &&
	   ( GenJetPt[corrJets(0)] > 0. && GenJetPt[corrJets(1)] > 0. ) &&
	   ( corrJets.pt(0)/GenJetPt[corrJets(0)] > 2.5 || corrJets.pt(1)/GenJetPt[corrJets(1)] > 2.5 ) ) {
	++nLargeResp;
	//std::cout << "WARNING: Large response, omitting event" << std::endl;
	continue;
      }


      // +++++ Data-based dijet selection ++++++++++++++++++++++++++++
      if( JetID[corrJets(0)] && JetID[corrJets(1)] ) {
	double deltaPhi = std::abs(TVector2::Phi_mpi_pi(JetPhi[corrJets(0)]-JetPhi[corrJets(1)]));
	bool passDeltaPhi = deltaPhi > minDeltaPhi_;
	double ptAve = 0.5*(corrJets.pt(0)+corrJets.pt(1));
	double pt2 = NObjJet > 2 ? corrJets.pt(2) : 0.;
	double ptrel = pt2/ptAve;
	double eta0 = JetEta[corrJets(0)];
	double eta1 = JetEta[corrJets(1)];
      
	if( ptrel <= binAdmin.ptSoftMax(binAdmin.nPtSoftBins()-1) ) {
	  unsigned int minPt3RelBin = 0;
	  while( ptrel > binAdmin.ptSoftMax(minPt3RelBin) ) ++minPt3RelBin;
      
	  unsigned int etaBin = 9999;
	  if( binAdmin.findSameEtaBin(eta0,eta1,etaBin) ) { // Same eta bin
	
	    unsigned int ptAveBin = 9999;
	    if( binAdmin.findPtBin(ptAve,etaBin,ptAveBin) ) { // PtAve bin
	  
	      if( !Parameters::runOnData || *(hltDecisions.at(etaBin).at(ptAveBin)) ) {
	    
		if( passDeltaPhi ) {
		  // Compute asymmetry
		  double ptAbsAsym = corrJets.pt(0)+corrJets.pt(1);
		  if( ptAbsAsym > 0 ) ptAbsAsym = (corrJets.pt(0) - corrJets.pt(1)) / ptAbsAsym;
		  double ptAsym = ( rand_->Uniform(0.,1.)>0.5 ? -1. : 1. ) * ptAbsAsym;
	      
		  // Loop over all pt3RelBin 
		  for(unsigned int pt3RelBin = minPt3RelBin; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
		
		    hPtAbsAsym.at(pt3RelBin).at(etaBin).at(ptAveBin)->Fill(ptAbsAsym,weight);
		    hPtAsym.at(pt3RelBin).at(etaBin).at(ptAveBin)->Fill(ptAsym,weight);
		    hPtAve.at(pt3RelBin).at(etaBin).at(ptAveBin)->Fill(ptAve,weight);
		  } // End of pt3Rel bin loop
	      
		} else {
		  ++nDeltaPhi;
		}
	    
	      } else {
		++nHlt;	    
	      }
	  
	    } else { // End of ptAve bin loop
	      ++nPtAve;
	    }
	
	  } else { // End of same eta bin loop
	    ++nEta;
	  }

	} else {
	  ++nJetID;
	}

      } else {
	++nPt3Rel;
      }



      // +++++ MC-truth dijet selection ++++++++++++++++++++++++++++++

      if( !Parameters::runOnData && NobjGenJet > 1 ) {

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
	double ptGenRel = ptGen2/ptGenAve;
	double pt0 = JetCorrL1[recoIdx0]*JetCorrL2L3[recoIdx0]*JetPt[recoIdx0];
	double pt1 = JetCorrL1[recoIdx1]*JetCorrL2L3[recoIdx1]*JetPt[recoIdx1];
	// Response of leading two gen-jets
	double r0 = 0.;
	double r1 = 0.;
	if( ptGen0 > 0. && ptGen1 > 0. ) {
	  r0 = pt0/ptGen0;
	  r1 = pt1/ptGen1;
	}
	// DeltaR between gen-jet and matched reco-jet
	double dr0 = util::deltaR(JetEta[recoIdx0],GenJetColEta[0],JetPhi[recoIdx0],GenJetColPhi[0]);
	double dr1 = util::deltaR(JetEta[recoIdx1],GenJetColEta[1],JetPhi[recoIdx1],GenJetColPhi[1]);
	bool closeMatch = ( dr0 < maxDeltaR_ && dr1 < maxDeltaR_ );
	// Do the matched reco-jets pass the jet id?
	bool passJetID = ( JetID[recoIdx0] && JetID[recoIdx1] );


	if( ptGenRel <= binAdmin.ptSoftMax(binAdmin.nPtSoftBins()-1) ) {
	  unsigned int minPt3RelBin = 0;
	  while( ptGenRel > binAdmin.ptSoftMax(minPt3RelBin) ) ++minPt3RelBin;
	
	  double deltaPhiGen = std::abs(TVector2::Phi_mpi_pi(GenJetColPhi[0]-GenJetColPhi[1]));
	  if( deltaPhiGen > minDeltaPhi_ && passJetID && closeMatch ) { 
	    
	    for(unsigned int pt3RelBin = minPt3RelBin; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {

	      unsigned int etaGenBin = 9999;
	      unsigned int ptGenBin = 9999;
	      if( binAdmin.findEtaBin(GenJetColEta[0],etaGenBin) ) {
		if( binAdmin.findPtBin(ptGen0,etaGenBin,ptGenBin) ) {
		  hResp.at(pt3RelBin).at(etaGenBin).at(ptGenBin)->Fill(r0,weight);
		}
	      }
	      if( binAdmin.findEtaBin(GenJetColEta[1],etaGenBin) ) {
		if( binAdmin.findPtBin(ptGen1,etaGenBin,ptGenBin) ) {
		  hResp.at(pt3RelBin).at(etaGenBin).at(ptGenBin)->Fill(r1,weight);
		}
	      }
	    } // End of loop over pt3rel bins

	  } // End if deltaPhi and jetID and DeltaR

	}	// End if valid pt3rel bin

      } // End if is MC and has >=2 genjets
      
    } // End of loop over entries
  
    delete chain;


    std::cout << "Writing output" << std::endl;

    TFile* resultFile = new TFile((Parameters::runOnData ? Parameters::histFileData : Parameters::histFileMC),"RECREATE");
  
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      TDirectory* dirEtaBin = resultFile->mkdir("Eta"+util::toTString(etaBin));
      
      for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
	TDirectory* dirPtBin = dirEtaBin->mkdir("Pt"+util::toTString(ptBin));

	for(unsigned int pt3RelBin = 0; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
	  TDirectory* outDir = dirPtBin->mkdir("PtSoft"+util::toTString(pt3RelBin));

	  outDir->WriteTObject(hPtAbsAsym.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hPtAsym.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hPtAve.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hResp.at(pt3RelBin).at(etaBin).at(ptBin));

	} // End of loop over pt3Rel bins
      } // End of loop over pt bins
    } // End of loop over eta bins

    resultFile->Close();
    delete resultFile;
  
  
  
    // ++++ Print status ++++++++++++++++++++++++++++++++++++++++++
    std::cout << "Done processing " << nEntries << " events from file list '" << inFileList << "'" << std::endl;
  
    std::cout << "Selected " << std::endl;
    std::cout << "  " << (nEntries -= nMaxNJet ) << " events with <= " << maxNJet << " jets " << std::endl;
    std::cout << "  " << (nEntries -= nDijets ) << " events with >= 2 jets " << std::endl;
    std::cout << "  " << (nEntries -= nJetID ) << " events with Jet(1,2) passing loose JetID cuts" << std::endl;
    if( !Parameters::runOnData ) std::cout << "  " << (nEntries -= nLargeResp ) << " events with good response" << std::endl;
    std::cout << "  " << (nEntries -= nPt3Rel ) << " events with nPt3Rel < " << binAdmin.ptSoftMax(binAdmin.nPtSoftBins()-1) << std::endl;
    std::cout << "  " << (nEntries -= nEta ) << " events with Eta(1,2) within binning" << std::endl;
    std::cout << "  " << (nEntries -= nPtAve ) << " events with PtAve within binning" << std::endl;
    if( Parameters::runOnData ) std::cout << "  " << (nEntries -= nHlt ) << " events passing HLT" << std::endl;
    std::cout << "  " << (nEntries -= nDeltaPhi ) << " events with |DeltaPhi(1,2)| > " << minDeltaPhi_ << std::endl;
  
  

    // ++++ Clean up +++++++++++++++++++++++++++++++++++++++++++++
  
    for(unsigned int pt3RelBin = 0; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
      for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
	for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
	  delete hPtAbsAsym.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hPtAsym.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hPtAve.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hResp.at(pt3RelBin).at(etaBin).at(ptBin);
	}
      }
    }

  }



  // --------------------------------------------------
  TChain *createTChain(const TString &fileListName) {
    std::cout << "  Getting trees from ntuples in '" << fileListName << "'" << std::endl;
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
      std::cerr << "  ERROR opening file '" << fileListName << "'\n";
      exit(1);
    }
    filelist.close();
    std::cout << "  Done. Added " << nOpenedFiles << " trees." << std::endl;

    return chain;
  }
}

