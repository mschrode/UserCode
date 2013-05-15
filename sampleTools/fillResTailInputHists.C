// $Id: fillResTailInputHists.C,v 1.1 2013/05/08 13:23:32 mschrode Exp $
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
#include "HistNames.h"
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
  void fillResTailInputHists(const Parameters &par, const sampleTools::BinningAdmin &binAdmin, bool runOnData, int nEvents);




  // ++++ FUNCTION IMPLEMENTATIONS +++++++++++++++++++++++++++++

  // --------------------------------------------------
  void fillResTailInputHists(const Parameters &par, const sampleTools::BinningAdmin &binAdmin, bool runOnData, int nEvents) {
  
    // ++++ Set parameters +++++++++++++++++++++++++++++++++++++++
  
    std::cout << "Initializing parameters" << std::endl;

    const TString inFileList = runOnData ? par.inFileListData() : par.inFileListMC();

    WeightProducer::init(par.qcdSampleType(),par.puScenario(),util::FileOps::readTH1(par.puHistFile(),"pileup"));

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

    std::cout << "Preparing histograms" << std::endl;

    HistNames hName;

    // Plots in analysis binning

    // Asymmetry distributions
    std::vector< std::vector< std::vector< TH1* > > > hPtAbsAsym;	// [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector< std::vector< TH1* > > > hPtAbsAsymGenBin;	// [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector< std::vector< TH1* > > > hPtAsym; // [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector< std::vector< TH1* > > > hPtAsymGenBin; // [pt3RelBin][etaBin][ptAveBin]

    // Control plots
    std::vector< std::vector< std::vector<TH1*> > > hPtAve; // [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector< std::vector<TH1*> > > hPtAveGenBin; // [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector< std::vector<TH1*> > > hRespGenBin; // [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector< std::vector<TH1*> > > hPt1; // [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector< std::vector<TH1*> > > hPt2; // [pt3RelBin][etaBin][ptAveBin]

    // N-1 plots of cut variabels
    std::vector< std::vector< std::vector<TH1*> > > hDeltaPhiNMin1; // [pt3RelBin][etaBin][ptAveBin]
    std::vector< std::vector<TH1*> > hPtAveNMin1;     // [pt3RelBin][etaBin]
    std::vector< std::vector<TH1*> > hPt3RelNMin1;    // [etaBin][ptAveBin]
  
    // Loop over pt3RelBins
    for(unsigned int pt3RelBin = 0; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {

      std::vector< std::vector< TH1* > > hPtAbsAsymTmp2D;
      std::vector< std::vector< TH1* > > hPtAbsAsymGenBinTmp2D;
      std::vector< std::vector< TH1* > > hPtAsymTmp2D;
      std::vector< std::vector< TH1* > > hPtAsymGenBinTmp2D;
      std::vector< std::vector< TH1* > > hPtAveTmp2D;
      std::vector< std::vector< TH1* > > hPtAveGenBinTmp2D;
      std::vector< std::vector< TH1* > > hRespGenBinTmp2D;
      std::vector< std::vector< TH1* > > hPt1Tmp2D;
      std::vector< std::vector< TH1* > > hPt2Tmp2D;
      std::vector< std::vector< TH1* > > hDeltaPhiNMin1Tmp2D;

      std::vector<TH1*> hPtAveNMin1Tmp; 
    

      // Loop over etaBins
      for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
	TString binId = 
	  "_Eta"+util::toTString(etaBin) +
	  "_PtSoft"+util::toTString(pt3RelBin);
	
	TString title = 
	  util::toTString(binAdmin.etaMin(etaBin)) + 
	  " < |#eta| < " + 
	  util::toTString(binAdmin.etaMax(etaBin)) +
	  ",  " + 
	  util::LabelFactory::pt3RelCut(binAdmin.ptSoftMax(pt3RelBin));

	TH1* h = new TH1D(hName.ptAve()+binId,title+";p^{ave}_{T} (GeV);Events",800,0,2000);
	h->Sumw2();
	hPtAveNMin1Tmp.push_back(h);


	std::vector< TH1* > hPtAbsAsymTmp;
	std::vector< TH1* > hPtAbsAsymGenBinTmp;
	std::vector< TH1* > hPtAsymTmp;
	std::vector< TH1* > hPtAsymGenBinTmp;
	std::vector< TH1* > hPtAveTmp;
	std::vector< TH1* > hPtAveGenBinTmp;
	std::vector< TH1* > hRespGenBinTmp;
	std::vector< TH1* > hPt1Tmp;
	std::vector< TH1* > hPt2Tmp;
	std::vector< TH1* > hDeltaPhiNMin1Tmp;

	// Loop over ptAveBins
	for(unsigned int ptAveBin = 0; ptAveBin < binAdmin.nPtBins(etaBin); ++ptAveBin) {

	  // Plots per (pt3Rel,eta,ptAve) bin
	  binId = 
	    "_Eta"+util::toTString(etaBin) +
	    "_Pt"+util::toTString(ptAveBin) +
	    "_PtSoft"+util::toTString(pt3RelBin);

	  title = 
	    util::toTString(binAdmin.etaMin(etaBin)) + 
	    " < |#eta| < " + 
	    util::toTString(binAdmin.etaMax(etaBin)) +
	    ",  " +
	    util::toTString(binAdmin.ptMin(etaBin,ptAveBin)) +
	    " < p^{ave}_{T} < " +
	    util::toTString(binAdmin.ptMax(etaBin,ptAveBin)) +	
	    " GeV,  " + util::LabelFactory::pt3RelCut(binAdmin.ptSoftMax(pt3RelBin));

	  TString genTitle = 
	    util::toTString(binAdmin.etaMin(etaBin)) + 
	    " < |#eta^{gen}| < " + 
	    util::toTString(binAdmin.etaMax(etaBin)) +
	    ",  " +
	    util::toTString(binAdmin.ptMin(etaBin,ptAveBin)) +
	    " < p^{gen,ave}_{T} < " +
	    util::toTString(binAdmin.ptMax(etaBin,ptAveBin)) +	
	    " GeV,  " + util::LabelFactory::pt3RelGenCut(binAdmin.ptSoftMax(pt3RelBin));

	  h = new TH1D(hName.ptAsymAbs()+binId,title+";|p_{T} asymmetry|;Events",150,-1.,1.);
	  h->Sumw2();
	  hPtAbsAsymTmp.push_back(h);

	  h = new TH1D(hName.ptAsymAbsGenBin()+binId,genTitle+";|p_{T} asymmetry|;Events",150,-1.,1.);
	  h->Sumw2();
	  hPtAbsAsymGenBinTmp.push_back(h);

	  h = new TH1D(hName.ptAsym()+binId,title+";p_{T} asymmetry;Events",101,-1.,1.);
	  h->Sumw2();
	  hPtAsymTmp.push_back(h);

	  h = new TH1D(hName.ptAsymGenBin()+binId,genTitle+";p_{T} asymmetry;Events",101,-1.,1.);
	  h->Sumw2();
	  hPtAsymGenBinTmp.push_back(h);

	  h = new TH1D(hName.ptAve()+binId,title+";p^{ave}_{T} (GeV);Events",800,0,2000);
	  h->Sumw2();
	  hPtAveTmp.push_back(h);

	  h = new TH1D(hName.ptAveGenBin()+binId,genTitle+";p^{ave}_{T} (GeV);Events",800,0,2000);
	  h->Sumw2();
	  hPtAveGenBinTmp.push_back(h);

	  h = new TH1D(hName.pt1()+binId,genTitle+";p_{T,1} (GeV);Events",800,0,2000);
	  h->Sumw2();
	  hPt1Tmp.push_back(h);

	  h = new TH1D(hName.pt2()+binId,genTitle+";p_{T,2} (GeV);Events",800,0,2000);
	  h->Sumw2();
	  hPt2Tmp.push_back(h);

	  h = new TH1D(hName.respGenBin()+binId,genTitle+";Response;Jets",120,0,2);
	  h->Sumw2();
	  hRespGenBinTmp.push_back(h);

	  h = new TH1D(hName.deltaPhi()+binId,title+";#Delta#phi;Events",320,0.,3.2);
	  h->Sumw2();
	  hDeltaPhiNMin1Tmp.push_back(h);
	}	// End of loop over ptAveBins

	hPtAbsAsymTmp2D.push_back(hPtAbsAsymTmp);
	hPtAbsAsymGenBinTmp2D.push_back(hPtAbsAsymGenBinTmp);
	hPtAsymTmp2D.push_back(hPtAsymTmp);
	hPtAsymGenBinTmp2D.push_back(hPtAsymGenBinTmp);
	hPtAveTmp2D.push_back(hPtAveTmp);
	hPtAveGenBinTmp2D.push_back(hPtAveGenBinTmp);
	hPt1Tmp2D.push_back(hPt1Tmp);
	hPt2Tmp2D.push_back(hPt2Tmp);
	hRespGenBinTmp2D.push_back(hRespGenBinTmp);
	hDeltaPhiNMin1Tmp2D.push_back(hDeltaPhiNMin1Tmp);

      } // End of loop over etaBin

      hPtAbsAsym.push_back(hPtAbsAsymTmp2D);
      hPtAbsAsymGenBin.push_back(hPtAbsAsymGenBinTmp2D);
      hPtAsym.push_back(hPtAsymTmp2D);
      hPtAsymGenBin.push_back(hPtAsymGenBinTmp2D);
      hPtAve.push_back(hPtAveTmp2D);
      hPtAveGenBin.push_back(hPtAveGenBinTmp2D);
      hPt1.push_back(hPt1Tmp2D);
      hPt2.push_back(hPt2Tmp2D);
      hRespGenBin.push_back(hRespGenBinTmp2D);
      hDeltaPhiNMin1.push_back(hDeltaPhiNMin1Tmp2D);

      hPtAveNMin1.push_back(hPtAveNMin1Tmp);

    } // End of loop over pt3RelBins


    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      std::vector<TH1*> hPt3RelNMin1Tmp;

      for(unsigned int ptAveBin = 0; ptAveBin < binAdmin.nPtBins(etaBin); ++ptAveBin) {
	TString binId = 
	  "_Eta"+util::toTString(etaBin) +
	  "_Pt"+util::toTString(ptAveBin);

	TString title = 
	  util::toTString(binAdmin.etaMin(etaBin)) + 
	  " < |#eta| < " + 
	  util::toTString(binAdmin.etaMax(etaBin)) +
	  ",  " +
	  util::toTString(binAdmin.ptMin(etaBin,ptAveBin)) +
	  " < p^{ave}_{T} < " +
	  util::toTString(binAdmin.ptMax(etaBin,ptAveBin));

	TH1* h = new TH1D(hName.pt3Rel()+binId,title+";"+util::LabelFactory::pt3Rel()+";Events",120,0.,1.2);
	h->Sumw2();
	hPt3RelNMin1Tmp.push_back(h);
      }

      hPt3RelNMin1.push_back(hPt3RelNMin1Tmp);
    }



    // Different binning
    // Helper histograms to access bins
    std::vector<double> ptBinEdges(26);
    util::HistOps::equidistLogBins(ptBinEdges,ptBinEdges.size()-1,4.,2700.);
    TH1* hPtGenBins = new TH1D("hPtGenBins",";p^{gen}_{T} (GeV)",ptBinEdges.size()-1,&(ptBinEdges.front()));

    // MC-truth response and control plots
    std::vector< std::vector<TH2*> > hRespVsPtGen; // [pt3RelBin][etaGenBin]
    std::vector< std::vector<TH2*> > hRespVsEtaGen; // [pt3RelBin][ptGenBin]

    std::vector<TH2*> hRespVsPtGenInclTmp;
    std::vector<TH2*> hRespVsEtaGenInclTmp;

    // Loop over pt3RelBins
    for(unsigned int pt3RelBin = 0; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {

      std::vector<TH2*> hRespVsPtGenTmp;
      for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
	TString binId = 
	  "_Eta"+util::toTString(etaBin) +
	  "_PtSoft"+util::toTString(pt3RelBin);
	
	TString title = 
	  util::toTString(binAdmin.etaMin(etaBin)) + 
	  " < |#eta^{gen}| < " + 
	  util::toTString(binAdmin.etaMax(etaBin)) +
	  ",  " + util::LabelFactory::pt3RelGenCut(binAdmin.ptSoftMax(pt3RelBin));

	TH2* h = new TH2D(hName.respVsPtGen()+binId,title+";p^{gen}_{T} (GeV);Response",
			  ptBinEdges.size()-1,&(ptBinEdges.front()),201,0.,2.);
	h->SetNdivisions(505);
	h->Sumw2();
	hRespVsPtGenTmp.push_back(h);

	if( pt3RelBin == 0 ) {
	  binId = "_Eta"+util::toTString(etaBin);
	  title = 
	    util::toTString(binAdmin.etaMin(etaBin)) + 
	    " < |#eta^{gen}| < " + 
	    util::toTString(binAdmin.etaMax(etaBin));
	  h = new TH2D(hName.respVsPtGen()+binId,title+";p^{gen}_{T} (GeV);Response",
		       ptBinEdges.size()-1,&(ptBinEdges.front()),201,0.,2.);
	  h->SetNdivisions(505);
	  h->Sumw2();
	  hRespVsPtGenInclTmp.push_back(h);
	}
      }
      hRespVsPtGen.push_back(hRespVsPtGenTmp);


      std::vector<TH2*> hRespVsEtaGenTmp;
      for(int ptGenBin = 0; ptGenBin < hPtGenBins->GetNbinsX(); ++ptGenBin) {
	TString binId = 
	  "_Pt"+util::toTString(ptGenBin) +
	  "_PtSoft"+util::toTString(pt3RelBin);
	
	TString title = 
	  util::toTString(hPtGenBins->GetXaxis()->GetBinLowEdge(ptGenBin+1)) +
	  " < p^{gen}_{T} < " +
	  util::toTString(hPtGenBins->GetXaxis()->GetBinUpEdge(ptGenBin+1)) +
	  ",  " + util::LabelFactory::pt3RelGenCut(binAdmin.ptSoftMax(pt3RelBin));

	TH2* h = new TH2D(hName.respVsEtaGen()+binId,title+";#eta^{gen};Response",
			  52,-5.2,5.2,201,0.,2.);
	h->SetNdivisions(505);
	h->Sumw2();
	hRespVsEtaGenTmp.push_back(h);

	if( pt3RelBin == 0 ) {
	  binId = "_Pt"+util::toTString(ptGenBin);	
	  title = 
	    util::toTString(hPtGenBins->GetXaxis()->GetBinLowEdge(ptGenBin+1)) +
	    " < p^{gen}_{T} < " +
	    util::toTString(hPtGenBins->GetXaxis()->GetBinUpEdge(ptGenBin+1));
	  h = new TH2D(hName.respVsEtaGen()+binId,title+";#eta^{gen};Response",
		       52,-5.2,5.2,201,0.,2.);
	  h->SetNdivisions(505);
	  h->Sumw2();
	  hRespVsEtaGenInclTmp.push_back(h);
	}
      }
      hRespVsEtaGen.push_back(hRespVsEtaGenTmp);
    }
    hRespVsPtGen.push_back(hRespVsPtGenInclTmp);
    hRespVsEtaGen.push_back(hRespVsEtaGenInclTmp);




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
    Float_t         GenJetPt[maxNJet];   //[NobjJet]
    Float_t         GenJetPhi[maxNJet];   //[NobjJet]
    Float_t         GenJetEta[maxNJet];   //[NobjJet]
    Int_t           NobjGenJet;
    Float_t         GenJetColPt[maxNJet];   //[NobjGenJet]
    Float_t         GenJetColPhi[maxNJet];   //[NobjGenJet]
    Float_t         GenJetColEta[maxNJet];   //[NobjGenJet]
    Int_t           GenJetColJetIdx[maxNJet];   //[NobjGenJet]
    Float_t         GenEvtScale;
    Int_t           PUMCNumVtx = 0;


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
    if( nEvents > 0 && nEvents <= nEntries ) nEntries = nEvents;

    for(int i = 0; i < nEntries; ++i) {
      if( i%50000 == 0 || i == nEntries-1 ) {
	std::cout << "  Processed " << i << " events" << std::endl;
      }

      chain->GetEntry(i);

      if( NObjJet > maxNJet ) {
	std::cerr << "WARNING: nObjJet = " << NObjJet << " > " << maxNJet << ". Skipping event!\n";
	++nMaxNJet;
	continue;
      }


      // Weights
      if( !runOnData ) {
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
      if ( !runOnData &&
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

	bool passesAnalysisSelection = false;
	unsigned int etaBin = 9999;
	unsigned int ptAveBin = 9999;
      
	unsigned int minPt3RelBin = 0;

	if( ptrel > binAdmin.ptSoftMax() ) {
	  minPt3RelBin = binAdmin.nPtSoftBins();
	} else {
	  while( ptrel > binAdmin.ptSoftMax(minPt3RelBin) ) ++minPt3RelBin;
	}
      
	if( binAdmin.findSameEtaBin(eta0,eta1,etaBin) ) { // Same eta bin
	  
	  if( binAdmin.findPtBin(ptAve,etaBin,ptAveBin) ) { // PtAve bin
	    
	    if( !runOnData || *(hltDecisions.at(etaBin).at(ptAveBin)) ) {
	      
	      // Analysis selection
	      if( passDeltaPhi ) {
		passesAnalysisSelection = true;

		if( minPt3RelBin < binAdmin.nPtSoftBins() ) {

		  // Compute asymmetry
		  double ptAbsAsym = corrJets.pt(0)+corrJets.pt(1);
		  if( ptAbsAsym > 0 ) ptAbsAsym = (corrJets.pt(0) - corrJets.pt(1)) / ptAbsAsym;
		  double ptAsym = ( rand_->Uniform(0.,1.)>0.5 ? -1. : 1. ) * ptAbsAsym;
		  
		  // Loop over all pt3RelBin 
		  for(unsigned int pt3RelBin = minPt3RelBin; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
		    hPtAbsAsym.at(pt3RelBin).at(etaBin).at(ptAveBin)->Fill(ptAbsAsym,weight);
		    hPtAsym.at(pt3RelBin).at(etaBin).at(ptAveBin)->Fill(ptAsym,weight);
		    hPtAve.at(pt3RelBin).at(etaBin).at(ptAveBin)->Fill(ptAve,weight);
		    
		    hPt1.at(pt3RelBin).at(etaBin).at(ptAveBin)->Fill(corrJets.pt(0),weight);
		    hPt2.at(pt3RelBin).at(etaBin).at(ptAveBin)->Fill(corrJets.pt(1),weight);

		    // N-1 plot for ptAve
		    hPtAveNMin1.at(pt3RelBin).at(etaBin)->Fill(ptAve,weight);

		  } // End of pt3Rel bin loop

		} else {
		  ++nPt3Rel;
		} // End of pt3rel bins selection
		
	      } else {
		++nDeltaPhi;
	      }
		
	      // N-1 plot
	      for(unsigned int pt3RelBin = minPt3RelBin; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
		hDeltaPhiNMin1.at(pt3RelBin).at(etaBin).at(ptAveBin)->Fill(deltaPhi,weight);
	      }		
	      
	    } else {
	      ++nHlt;	    
	    }	// End of HLT selection
	  
	  } else {
	    ++nPtAve;
	  } // End of ptAve bin selection

	} else {
	  ++nEta;
	} // End of if jets in same eta bin switch

	// N-1 plot for pt3Rel
	if( passesAnalysisSelection ) {
	  hPt3RelNMin1.at(etaBin).at(ptAveBin)->Fill(ptrel,weight);
	}

      } else {
	++nJetID;
      }	// End of jetID selection


      // +++++ MC-truth dijet selection ++++++++++++++++++++++++++++++

      if( !runOnData && NobjGenJet > 1 ) {

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
	
	  if( passJetID && closeMatch ) { 
	    double deltaPhiGen = std::abs(TVector2::Phi_mpi_pi(GenJetColPhi[0]-GenJetColPhi[1]));
	    unsigned int etaGenBin = 9999;
	    unsigned int ptGenBin = 9999;

	    for(unsigned int pt3RelBin = minPt3RelBin; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
	      // Analysis cuts
	      if( deltaPhiGen > minDeltaPhi_ ) {
		
		// Asymmetry in gen-bins
		if( binAdmin.findSameEtaBin(GenJetColEta[0],GenJetColEta[1],etaGenBin) ) { // Same eta bin
		  if( binAdmin.findPtBin(ptGenAve,etaGenBin,ptGenBin) ) { // PtGenAve bin
		    // Compute asymmetry
		    double ptAbsAsym = pt0+pt1;
		    if( ptAbsAsym > 0 ) ptAbsAsym = (pt0-pt1) / ptAbsAsym;
		    double ptAsym = ( rand_->Uniform(0.,1.)>0.5 ? -1. : 1. ) * ptAbsAsym;
		    hPtAbsAsymGenBin.at(pt3RelBin).at(etaGenBin).at(ptGenBin)->Fill(ptAbsAsym,weight);
		    hPtAsymGenBin.at(pt3RelBin).at(etaGenBin).at(ptGenBin)->Fill(ptAsym,weight);
		    hPtAveGenBin.at(pt3RelBin).at(etaGenBin).at(ptGenBin)->Fill(ptGenAve,weight);
		  }
		}
		
		// Response and MC truth; needs only one jet in bin
		if( binAdmin.findEtaBin(GenJetColEta[0],etaGenBin) ) {
		  if( binAdmin.findPtBin(ptGen0,etaGenBin,ptGenBin) ) {
		    hRespGenBin.at(pt3RelBin).at(etaGenBin).at(ptGenBin)->Fill(r0,weight);
		  }
		}
		if( binAdmin.findEtaBin(GenJetColEta[1],etaGenBin) ) {
		  if( binAdmin.findPtBin(ptGen1,etaGenBin,ptGenBin) ) {
		    hRespGenBin.at(pt3RelBin).at(etaGenBin).at(ptGenBin)->Fill(r1,weight);
		  }
		}

	      } // End if deltaPhi


	      // MC-truth response histograms
	      if( binAdmin.findEtaBin(GenJetColEta[0],etaGenBin) ) {
		hRespVsPtGen.at(pt3RelBin).at(etaGenBin)->Fill(ptGen0,r0,wPU);
	      }
	      if( binAdmin.findEtaBin(GenJetColEta[1],etaGenBin) ) {
		hRespVsPtGen.at(pt3RelBin).at(etaGenBin)->Fill(ptGen1,r1,wPU);
	      }

	      ptGenBin = hPtGenBins->FindBin(ptGen0)-1;
	      if( ptGenBin < hRespVsEtaGen.at(pt3RelBin).size() ) {
		hRespVsEtaGen.at(pt3RelBin).at(ptGenBin)->Fill(GenJetColEta[0],r0,wPU);
	      }
	      ptGenBin = hPtGenBins->FindBin(ptGen1)-1;
	      if( ptGenBin < hRespVsEtaGen.at(pt3RelBin).size() ) {
		hRespVsEtaGen.at(pt3RelBin).at(ptGenBin)->Fill(GenJetColEta[1],r1,wPU);
	      }
	    } // End of loop over pt3rel bins

	    // MC-truth response histograms
	    // Incls in pt3
	    unsigned int pt3RelBin = binAdmin.nPtSoftBins();
	    if( binAdmin.findEtaBin(GenJetColEta[0],etaGenBin) ) {
	      hRespVsPtGen.at(pt3RelBin).at(etaGenBin)->Fill(ptGen0,r0,wPU);
	    }
	    if( binAdmin.findEtaBin(GenJetColEta[1],etaGenBin) ) {
	      hRespVsPtGen.at(pt3RelBin).at(etaGenBin)->Fill(ptGen1,r1,wPU);
	    }
	    
	    ptGenBin = hPtGenBins->FindBin(ptGen0)-1;
	    if( ptGenBin < hRespVsEtaGen.at(pt3RelBin).size() ) {
	      hRespVsEtaGen.at(pt3RelBin).at(ptGenBin)->Fill(GenJetColEta[0],r0,wPU);
	    }
	    ptGenBin = hPtGenBins->FindBin(ptGen1)-1;
	    if( ptGenBin < hRespVsEtaGen.at(pt3RelBin).size() ) {
	      hRespVsEtaGen.at(pt3RelBin).at(ptGenBin)->Fill(GenJetColEta[1],r1,wPU);
	    }

	  } // End if jetID and DeltaR

	} // End if valid pt3rel bin

      } // End if is MC and has >=2 genjets

    } // End of loop over entries
  
    delete chain;


    std::cout << "Writing output" << std::endl;

    TFile* resultFile = new TFile((runOnData ? par.histFileData() : par.histFileMC()),"RECREATE");
  
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      TDirectory* dirEtaBin = resultFile->mkdir("Eta"+util::toTString(etaBin));
      
      for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
	TDirectory* dirPtBin = dirEtaBin->mkdir("Pt"+util::toTString(ptBin));

	for(unsigned int pt3RelBin = 0; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
	  TDirectory* outDir = dirPtBin->mkdir("PtSoft"+util::toTString(pt3RelBin));

	  outDir->WriteTObject(hPtAbsAsym.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hPtAbsAsymGenBin.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hPtAsym.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hPtAsymGenBin.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hPtAve.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hPtAveGenBin.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hRespGenBin.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hPt1.at(pt3RelBin).at(etaBin).at(ptBin));
	  outDir->WriteTObject(hPt2.at(pt3RelBin).at(etaBin).at(ptBin));
	} // End of loop over pt3Rel bins

      } // End of loop over pt bins

    } // End of loop over eta bins

    TDirectory* outDir = resultFile->mkdir("NMinus1Plots");
    for(unsigned int pt3RelBin = 0; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
      for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
	outDir->WriteTObject(hPtAveNMin1.at(pt3RelBin).at(etaBin));
	for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
	  outDir->WriteTObject(hDeltaPhiNMin1.at(pt3RelBin).at(etaBin).at(ptBin));
	}
      }
    }
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
	outDir->WriteTObject(hPt3RelNMin1.at(etaBin).at(ptBin));
      }
    }

    if( !runOnData ) {
      outDir = resultFile->mkdir("MCTruthResolution");
      for(unsigned int pt3RelBin = 0; pt3RelBin <= binAdmin.nPtSoftBins(); ++pt3RelBin) {
	for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
	  outDir->WriteTObject(hRespVsPtGen.at(pt3RelBin).at(etaBin));
	}
	for(unsigned int ptGenBin = 0; ptGenBin < hRespVsEtaGen.at(pt3RelBin).size(); ++ptGenBin) {
	  outDir->WriteTObject(hRespVsEtaGen.at(pt3RelBin).at(ptGenBin));
	}
      }
    }

    resultFile->Close();
    delete resultFile;
  
  
  
    // ++++ Print status ++++++++++++++++++++++++++++++++++++++++++
    std::cout << "Done processing " << nEntries << " events from file list '" << inFileList << "'" << std::endl;
  
    std::cout << "Selected " << std::endl;
    std::cout << "  " << (nEntries -= nMaxNJet ) << " events with <= " << maxNJet << " jets " << std::endl;
    std::cout << "  " << (nEntries -= nDijets ) << " events with >= 2 jets " << std::endl;
    if( !runOnData ) std::cout << "  " << (nEntries -= nLargeResp ) << " events with good response" << std::endl;
    std::cout << "  " << (nEntries -= nJetID ) << " events with Jet(1,2) passing loose JetID cuts" << std::endl;
    std::cout << "  " << (nEntries -= nEta ) << " events with Eta(1,2) within binning" << std::endl;
    std::cout << "  " << (nEntries -= nPtAve ) << " events with PtAve within binning" << std::endl;
    if( runOnData ) std::cout << "  " << (nEntries -= nHlt ) << " events passing HLT" << std::endl;
    std::cout << "  " << (nEntries -= nDeltaPhi ) << " events with |DeltaPhi(1,2)| > " << minDeltaPhi_ << std::endl;
    std::cout << "  " << (nEntries -= nPt3Rel ) << " events with nPt3Rel < " << binAdmin.ptSoftMax(binAdmin.nPtSoftBins()-1) << std::endl;
  
  

    // ++++ Clean up +++++++++++++++++++++++++++++++++++++++++++++
  
    for(unsigned int pt3RelBin = 0; pt3RelBin < binAdmin.nPtSoftBins(); ++pt3RelBin) {
      for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
	for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
	  delete hPtAbsAsym.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hPtAbsAsymGenBin.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hPtAsym.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hPtAsymGenBin.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hPtAve.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hPtAveGenBin.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hRespGenBin.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hPt1.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hPt2.at(pt3RelBin).at(etaBin).at(ptBin);
	  delete hDeltaPhiNMin1.at(pt3RelBin).at(etaBin).at(ptBin);
	}
	delete hPtAveNMin1.at(pt3RelBin).at(etaBin);
      }
    }
    for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
      for(unsigned int ptBin = 0; ptBin < binAdmin.nPtBins(etaBin); ++ptBin) {
	delete hPt3RelNMin1.at(etaBin).at(ptBin);
      }
    }
    for(unsigned int pt3RelBin = 0; pt3RelBin <= binAdmin.nPtSoftBins(); ++pt3RelBin) {
      for(unsigned int etaBin = 0; etaBin < binAdmin.nEtaBins(); ++etaBin) {
	delete hRespVsPtGen.at(pt3RelBin).at(etaBin);
      }
      for(unsigned int ptGenBin = 0; ptGenBin < hRespVsEtaGen.at(pt3RelBin).size(); ++ptGenBin) {
	delete hRespVsEtaGen.at(pt3RelBin).at(ptGenBin);
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

