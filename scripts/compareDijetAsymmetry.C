#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TString.h"
#include "TVector2.h"

#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/utils.h"
#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/HistOps.h"
#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/LabelFactory.h"
#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/StyleSettings.h"


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


void compareDijetAsymmetry(int nMaxEvts = -1, int prescale = 1) {
  util::StyleSettings::presentationNoTitle();

  const TString prefix = "JetMET_Run2010A-PromptReco-v4_614nb_Pt3RelCuts_Eta13-30_";
  const int maxNJet = 50;
  const double minEMF = 0.01;
  const int minJetN90Hits = 2;
  const double maxJetFHPD = 0.98;
  const double maxJetFRBX = 0.98;
  const double minDeltaPhi = 2.7;
  const double minEta = 1.3;
  const double maxEta = 3.;
  const double minPtAve = 20.;
  const double maxPt3Rel = 0.1;
  const double lumi = 0.614;
  const double lumiFrac = lumi/0.1;
  const int mcColor = 5;
  TRandom *rand = new TRandom(0);

  // 0: Data
  // 1: MC
  unsigned int nChains = 2;

  TChain *data = createTChain("/afs/naf.desy.de/user/m/mschrode/Kalibri/input/Kalibri_JetMET_Run2010A-PromptReco-v4_DCSONLY"); //614 nb
  TChain *mc = createTChain("/afs/naf.desy.de/user/m/mschrode/Kalibri/input/Kalibri_Spring10QCDDiJet_AK5Calo");


  std::vector<int> nTotal(nChains,0);
  std::vector<int> nPtAve(nChains,0);
  std::vector<int> nEta(nChains,0);
  std::vector<int> nDeltaPhi(nChains,0);
  std::vector<int> nPt3Rel(nChains,0);
  std::vector<int> nJetID(nChains,0);


  // Histograms
  std::vector<double> ptAveBinEdges;
  ptAveBinEdges.push_back(10.);
  ptAveBinEdges.push_back(20.);
  ptAveBinEdges.push_back(30.);
  ptAveBinEdges.push_back(40.);
  ptAveBinEdges.push_back(60.);
  ptAveBinEdges.push_back(80.);
  ptAveBinEdges.push_back(100.);
  ptAveBinEdges.push_back(120.);
  ptAveBinEdges.push_back(140.);
  ptAveBinEdges.push_back(170.);
  ptAveBinEdges.push_back(200.);
  ptAveBinEdges.push_back(250.);
  ptAveBinEdges.push_back(300.);
  ptAveBinEdges.push_back(350.);
  ptAveBinEdges.push_back(400.);
  ptAveBinEdges.push_back(1000.);

  std::vector<double> ptAveBinEdgesLog(50);
  util::HistOps::equidistLogBins(ptAveBinEdgesLog,ptAveBinEdgesLog.size()-1,ptAveBinEdges.front(),ptAveBinEdges.back());
  std::vector<double> ptJet3BinEdgesLog(50);
  util::HistOps::equidistLogBins(ptJet3BinEdgesLog,ptJet3BinEdgesLog.size()-1,1.,300.);

  std::vector< std::vector<TH1*> > hPtUncorr(nChains);
  std::vector< std::vector<TH1*> > hPtCorr(nChains);
  std::vector<TH2*> hPtJet1vs2Corr(nChains);
  std::vector<TH2*> hPtJet1vs2Uncorr(nChains);
  std::vector<TH2*> hPtJet1vs2CorrLog(nChains);
  std::vector<TH2*> hPtJet1vs2UncorrLog(nChains);
  std::vector<TH1*> hPtAveCorr(nChains);
  std::vector<TH1*> hPtAveCorrLog(nChains);
  std::vector<TH2*> hPtAsymVsPtAveUncorr(nChains);
  std::vector<TH2*> hPtAsymVsPtAveCorr(nChains);
  std::vector<TH2*> hPtBiasAsymVsPtAveCorr(nChains);
  std::vector<TH2*> hPtAsymVsPtJet3Uncorr(nChains);
  std::vector<TH1*> hDeltaPhi(nChains);
  std::vector<TH1*> hEta(nChains);
  std::vector<TH2*> hPtAsymVsEMF(nChains);
  for(size_t c = 0; c < nChains; ++c) {
    TH1 *h = 0;
    for(int j = 0; j < 3; ++j) {
      h = util::HistOps::createTH1D("hPtUncorr"+util::toTString(c)+"_"+util::toTString(j+1),50,ptAveBinEdges.front(),ptAveBinEdges.back(),"p_{T,"+util::toTString(j+1)+"}","GeV","jets",false);
      h->Sumw2();
      if( c == 0 ) h->SetMarkerStyle(20);
      else h->SetFillColor(mcColor);
      hPtUncorr[c].push_back(h);

      h = util::HistOps::createTH1D("hPtCorr"+util::toTString(c)+"_"+util::toTString(j+1),50,ptAveBinEdges.front(),ptAveBinEdges.back(),"Corrected p_{T,"+util::toTString(j+1)+"}","GeV","jets",false);
      h->Sumw2();
      if( c == 0 ) h->SetMarkerStyle(20);
      else h->SetFillColor(mcColor);
      hPtCorr[c].push_back(h);
    }
    h = util::HistOps::createTH1D("hPtAveCorr"+util::toTString(c),50,ptAveBinEdges.front(),ptAveBinEdges.back(),"Corrected p^{ave}_{T}","GeV","jets",false);
    h->Sumw2();
    if( c == 0 ) h->SetMarkerStyle(20);
    else h->SetFillColor(mcColor);
    hPtAveCorr[c] = h;

    h = util::HistOps::createTH1D("hPtAveCorrLog"+util::toTString(c),ptAveBinEdgesLog.size()-1,&(ptAveBinEdgesLog.front()),"Corrected p^{ave}_{T}","GeV","jets",false);
    h->Sumw2();
    if( c == 0 ) h->SetMarkerStyle(20);
    else h->SetFillColor(mcColor);
    hPtAveCorrLog[c] = h;

    h = util::HistOps::createTH1D("hDeltaPhi"+util::toTString(c),50,0.,M_PI,"|#Delta#phi_{12}|","","events",false);
    h->Sumw2();
    if( c == 0 ) h->SetMarkerStyle(20);
    else h->SetFillColor(mcColor);
    hDeltaPhi[c] = h;

    h = util::HistOps::createTH1D("hEta"+util::toTString(c),51,-5.1,5.1,"#eta_{1,2}","","jets",false);
    h->Sumw2();
    if( c == 0 ) h->SetMarkerStyle(20);
    else h->SetFillColor(mcColor);
    hEta[c] = h;

    TH2 *h2 = new TH2D("hPtAsymVsPtAveUncorr"+util::toTString(c),"",ptAveBinEdges.size()-1,&(ptAveBinEdges.front()),31,-1.,1.);
    h2->Sumw2();
    hPtAsymVsPtAveUncorr[c] = h2;

    h2 = new TH2D("hPtAsymVsPtAveCorr"+util::toTString(c),"",ptAveBinEdges.size()-1,&(ptAveBinEdges.front()),31,-1.,1.);
    h2->Sumw2();
    hPtAsymVsPtAveCorr[c] = h2;

    h2 = new TH2D("hPtBiasAsymVsPtAveCorr"+util::toTString(c),"",ptAveBinEdges.size()-1,&(ptAveBinEdges.front()),15,0.,1.);
    h2->Sumw2();
    hPtBiasAsymVsPtAveCorr[c] = h2;

    h2 = new TH2D("hPtAsymVsPtJet3Uncorr"+util::toTString(c),"p_{T,3} (GeV);Asymmetry",ptJet3BinEdgesLog.size()-1,&(ptJet3BinEdgesLog.front()),31,-1.,1.);
    h2->Sumw2();
    hPtAsymVsPtJet3Uncorr[c] = h2;

    hPtJet1vs2Uncorr[c] = new TH2D("hPtJet1vs2Uncorr"+util::toTString(c),";p_{T,1} (GeV);p_{T,2} (GeV)",50,0.,200.,50,0.,200.);

    hPtJet1vs2Corr[c] = new TH2D("hPtJet1vs2Corr"+util::toTString(c),";Corrected p_{T,1} (GeV);Corrected p_{T,2} (GeV)",50,0.,200.,50,0.,200.);

    hPtJet1vs2UncorrLog[c] = new TH2D("hPtJet1vs2UncorrLog"+util::toTString(c),";p_{T,1} (GeV);p_{T,2} (GeV)",ptAveBinEdgesLog.size()-1,&(ptAveBinEdgesLog.front()),ptAveBinEdgesLog.size()-1,&(ptAveBinEdgesLog.front()));

    hPtJet1vs2CorrLog[c] = new TH2D("hPtJet1vs2CorrLog"+util::toTString(c),";Corrected p_{T,1} (GeV);Corrected p_{T,2} (GeV)",ptAveBinEdgesLog.size()-1,&(ptAveBinEdgesLog.front()),ptAveBinEdgesLog.size()-1,&(ptAveBinEdgesLog.front()));

    hPtAsymVsEMF[c] = new TH2D("hPtAsymVsEMF"+util::toTString(c),";f_{em};Asymmetry",50,-1.5,1.5,31,-1.,1.);
  }


  // Read variables
  std::vector<JetIndex*> jIdx;
  jIdx.resize(8,0);
  for(unsigned int c = 0; c < nChains; ++c) {
    std::cout << "\nReading events from chain " << c << std::endl;
    TChain *chain = 0;
    if     ( c == 0 ) chain = data;
    else if( c == 1 ) chain = mc;

    if( chain ) {
      // Init read quantities
      float weight = 1.;
      int nObjJet = 0;
      float jetPt[maxNJet];
      float jetEta[maxNJet];
      float jetPhi[maxNJet];
      float jetEMF[maxNJet];
      int jetN90Hits[maxNJet];
      float jetFHPD[maxNJet];
      float jetFRBX[maxNJet];
      float jetCorrL2L3[maxNJet];

      // Set branch addresses
      if( c > 0 ) chain->SetBranchAddress("Weight",&weight);
      chain->SetBranchAddress("NobjJet",&nObjJet);
      chain->SetBranchAddress("JetPt",jetPt);
      chain->SetBranchAddress("JetEta",jetEta);
      chain->SetBranchAddress("JetPhi",jetPhi);
      chain->SetBranchAddress("JetEMF",jetEMF);
      chain->SetBranchAddress("JetN90Hits",jetN90Hits);
      chain->SetBranchAddress("JetFHPD",jetFHPD);
      chain->SetBranchAddress("JetFRBX",jetFRBX);
      chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  
      // Loop over tree entries and fill histograms
      int nEntries = chain->GetEntries();
      if( nMaxEvts > 0 && nEntries > nMaxEvts ) nEntries = nMaxEvts;
      for(int n = 0; n < nEntries; n = n + prescale) {
	if( n%100000 == 0 ) std::cout << " Entry " << n << std::endl;
	chain->GetEntry(n);

	if( nObjJet > maxNJet ) {
	  std::cerr << "WARNING: nObjJet = " << nObjJet << " > " << maxNJet << ". Skipping event.\n";
	  continue;
	}


	if( nObjJet < 2 ) continue;
	nTotal[c]++;

	if( c == 1 ) weight *= lumiFrac;

	// Uncorrected
	bool isGoodEvt = true;
	double ptAveCorr = 0.5*(jetCorrL2L3[0]*jetPt[0]+jetCorrL2L3[1]*jetPt[1]);
	if( 0.5*(jetPt[0]+jetPt[1]) < minPtAve ) isGoodEvt = false;
	else if( std::abs(jetEta[0]) < minEta || std::abs(jetEta[1]) < minEta ||
		 std::abs(jetEta[0]) > maxEta || std::abs(jetEta[1]) > maxEta ) isGoodEvt = false;
	else if( std::abs(TVector2::Phi_mpi_pi(jetPhi[0]-jetPhi[1])) < minDeltaPhi ) isGoodEvt = false;
	else if( (nObjJet>2 ? jetPt[2] : 0.) > maxPt3Rel*ptAveCorr ) isGoodEvt = false;
	else if( jetN90Hits[0] < minJetN90Hits || jetN90Hits[1] < minJetN90Hits ) isGoodEvt = false;
	else if( jetFHPD[0] > maxJetFHPD || jetFHPD[1] > maxJetFHPD ) isGoodEvt = false;
	else if( (std::abs(jetEta[0]) < 2.55 && jetEMF[0] < minEMF) ||
		 (std::abs(jetEta[1]) < 2.55 && jetEMF[1] < minEMF) ) isGoodEvt = false;
	else if( (std::abs(jetEta[0]) > 2.55 && jetEMF[0] < -0.9) ||
		 (std::abs(jetEta[1]) > 2.55 && jetEMF[1] < -0.9) ) isGoodEvt = false;
	else if( (std::abs(jetEta[0]) > 2.55 && jetPt[0] > 80. && jetEMF[0] > 1. ) ||
		 (std::abs(jetEta[1]) > 2.55 && jetPt[1] > 80. && jetEMF[1] > 1. ) ) isGoodEvt = false;

	if( isGoodEvt ) {
	  hPtUncorr[c][0]->Fill(jetPt[0],weight);
	  hPtUncorr[c][1]->Fill(jetPt[1],weight);
	  if( nObjJet > 2 ) hPtUncorr[c][2]->Fill(jetPt[2],weight);

	  double ptAsym = jetPt[0]+jetPt[1];
	  if( ptAsym > 0. ) {
	    ptAsym = (jetPt[0]-jetPt[1])/ptAsym;
	    if( rand->Uniform()>0.5 ) {
	      ptAsym *= -1.;
	      hPtJet1vs2Uncorr[c]->Fill(jetPt[1],jetPt[0],weight);
	      hPtJet1vs2UncorrLog[c]->Fill(jetPt[1],jetPt[0],weight);
	    } else { 
	      hPtJet1vs2Uncorr[c]->Fill(jetPt[0],jetPt[1],weight);
	      hPtJet1vs2UncorrLog[c]->Fill(jetPt[0],jetPt[1],weight);
	    }
	    hPtAsymVsPtAveUncorr[c]->Fill(0.5*(jetPt[0]+jetPt[1]),ptAsym,weight);
	    if( nObjJet > 2 ) hPtAsymVsPtJet3Uncorr[c]->Fill(jetPt[2],ptAsym,weight);
	    hPtAsymVsEMF[c]->Fill(jetEMF[0],ptAsym,weight);
	    hPtAsymVsEMF[c]->Fill(jetEMF[1],ptAsym,weight);
	  }
	}


 	// Corrected
 	isGoodEvt = true;

 	// Sort by corrected pt
 	unsigned int nJets = nObjJet;
 	if( jIdx.size() < nJets ) nJets = jIdx.size();
 	for(size_t i = 0; i < nJets; ++i) {
 	  jIdx[i] = new JetIndex(i,jetPt[i]*jetCorrL2L3[i]);
 	}
 	//std::sort(jIdx.begin(),jIdx.begin()+nJets,JetIndex::ptGreaterThan);

  	// Select events
	ptAveCorr = 0.5*(jIdx[0]->pt_+jIdx[1]->pt_);
  	if( ptAveCorr < minPtAve ) {
	  nPtAve[c]++;
	  isGoodEvt = false;
	}
  	else if( std::abs(jetEta[jIdx[0]->idx_]) < minEta || std::abs(jetEta[jIdx[1]->idx_]) < minEta ||
		 std::abs(jetEta[jIdx[0]->idx_]) > maxEta || std::abs(jetEta[jIdx[1]->idx_]) > maxEta ) {
	  nEta[c]++;
	  isGoodEvt = false;
	}
  	else if( std::abs(TVector2::Phi_mpi_pi(jetPhi[jIdx[0]->idx_]-jetPhi[jIdx[1]->idx_])) < minDeltaPhi ) {
	  nDeltaPhi[c]++;
	  isGoodEvt = false;
	}
	else if( (nObjJet>2 ? jetPt[2] : 0.) > maxPt3Rel*ptAveCorr ) {
	  nPt3Rel[c]++;
	  isGoodEvt = false;
	}
	else if( jetN90Hits[jIdx[0]->idx_] < minJetN90Hits || jetN90Hits[jIdx[1]->idx_] < minJetN90Hits ) {
	  nJetID[c]++;
	  isGoodEvt = false;
	}
	else if( jetFHPD[jIdx[0]->idx_] > maxJetFHPD || jetFHPD[jIdx[1]->idx_] > maxJetFHPD ) {
	  nJetID[c]++;
	  isGoodEvt = false;
	}
	else if( (std::abs(jetEta[jIdx[0]->idx_]) < 2.55 && jetEMF[jIdx[0]->idx_] < minEMF) ||
		 (std::abs(jetEta[jIdx[1]->idx_]) < 2.55 && jetEMF[jIdx[1]->idx_] < minEMF) ) {
	  nJetID[c]++;
	  isGoodEvt = false;
	}
	else if( (std::abs(jetEta[jIdx[0]->idx_]) > 2.55 && jetEMF[jIdx[0]->idx_] < -0.9) ||
		 (std::abs(jetEta[jIdx[1]->idx_]) > 2.55 && jetEMF[jIdx[1]->idx_] < -0.9) ) {
	  nJetID[c]++;
	  isGoodEvt = false;
	}
	else if( (std::abs(jetEta[jIdx[0]->idx_]) > 2.55 && jetPt[jIdx[0]->idx_] > 80. && jetEMF[jIdx[0]->idx_] > 1. ) ||
		 (std::abs(jetEta[jIdx[1]->idx_]) > 2.55 && jetPt[jIdx[1]->idx_] > 80. && jetEMF[jIdx[1]->idx_] > 1. ) ) {
	  nJetID[c]++;
	  isGoodEvt = false;
	}

  	if( isGoodEvt ) {
  	  hPtCorr[c][0]->Fill(jIdx[0]->pt_,weight);
  	  hPtCorr[c][1]->Fill(jIdx[1]->pt_,weight);
  	  if( nObjJet > 2 ) hPtCorr[c][2]->Fill(jIdx[2]->pt_,weight);
	  hPtAveCorr[c]->Fill(ptAveCorr,weight);
	  hPtAveCorrLog[c]->Fill(ptAveCorr,weight);
	  hDeltaPhi[c]->Fill(std::abs(TVector2::Phi_mpi_pi(jetPhi[jIdx[0]->idx_]-jetPhi[jIdx[1]->idx_])),weight);
	  hEta[c]->Fill(jetEta[jIdx[0]->idx_],weight);
	  hEta[c]->Fill(jetEta[jIdx[1]->idx_],weight);

  	  double ptAsym = jIdx[0]->pt_+jIdx[1]->pt_;
  	  if( ptAsym > 0. ) {
  	    ptAsym = (jIdx[0]->pt_-jIdx[1]->pt_)/ptAsym;
	    hPtBiasAsymVsPtAveCorr[c]->Fill(ptAveCorr,ptAsym,weight);
	    if( rand->Uniform()>0.5 ) {
	      ptAsym *= -1.;
	      hPtJet1vs2Corr[c]->Fill(jIdx[1]->pt_,jIdx[0]->pt_,weight);
	      hPtJet1vs2CorrLog[c]->Fill(jIdx[1]->pt_,jIdx[0]->pt_,weight);
	    } else {
	      hPtJet1vs2Corr[c]->Fill(jIdx[0]->pt_,jIdx[1]->pt_,weight);
	      hPtJet1vs2CorrLog[c]->Fill(jIdx[0]->pt_,jIdx[1]->pt_,weight);
	    }
  	    hPtAsymVsPtAveCorr[c]->Fill(ptAveCorr,ptAsym,weight);
  	  }
  	}
 	for(size_t i = 0; i < nJets; ++i) {
 	  delete jIdx[i];
 	}
      } // End of loop over entries
    }
  } // End of loop over chains


  // Get pt asymmetry per bin
  std::vector< std::vector<TH1*> > hPtAsymUncorr(nChains);
  std::vector< std::vector<TH1*> > hPtAsymCorr(nChains);
  std::vector< std::vector<TH1*> > hPtBiasAsymCorr(nChains);
  std::vector<TH1*> hPtAsymCorrWidthGauss(nChains);
  std::vector<TH1*> hPtAsymCorrWidthStd(nChains);
  for(unsigned int c = 0; c < nChains; ++c) {
    util::HistOps::fillSlices(hPtAsymVsPtAveUncorr[c],hPtAsymUncorr[c],"hPtAsymUncorr"+util::toTString(c));
    util::HistOps::fillSlices(hPtAsymVsPtAveCorr[c],hPtAsymCorr[c],"hPtAsymCorr"+util::toTString(c));
    util::HistOps::fillSlices(hPtBiasAsymVsPtAveCorr[c],hPtBiasAsymCorr[c],"hPtBiasAsymCorr"+util::toTString(c));

    hPtAsymCorrWidthGauss[c] = util::HistOps::createTH1D("hPtAsymCorrWidthGauss"+util::toTString(c),ptAveBinEdges.size()-1,&(ptAveBinEdges.front()),"p^{ave}_{T}","GeV","#sigma(Asymmetry)");
    if( c == 0 ) hPtAsymCorrWidthGauss[c]->SetMarkerStyle(20);
    else hPtAsymCorrWidthGauss[c]->SetFillColor(mcColor);
    
    hPtAsymCorrWidthStd[c] = util::HistOps::createTH1D("hPtAsymCorrWidthStd"+util::toTString(c),ptAveBinEdges.size()-1,&(ptAveBinEdges.front()),"p^{ave}_{T}","GeV","StdDev(Asymmetry)");
    if( c == 0 ) hPtAsymCorrWidthStd[c]->SetMarkerStyle(20);
    else hPtAsymCorrWidthStd[c]->SetFillColor(mcColor);

    for(size_t p = 0; p < hPtAsymUncorr[c].size(); ++p) {
      util::HistOps::setAxisTitles(hPtBiasAsymCorr[c][p],"Biased Asymmetry","","events",false);
      if( c == 0 ) hPtBiasAsymCorr[c][p]->SetMarkerStyle(20);
      else hPtBiasAsymCorr[c][p]->SetFillColor(mcColor);

      util::HistOps::setAxisTitles(hPtAsymUncorr[c][p],"Asymmetry","","events",false);
      if( c == 0 ) hPtAsymUncorr[c][p]->SetMarkerStyle(20);
      else hPtAsymUncorr[c][p]->SetFillColor(mcColor);
      
      util::HistOps::setAxisTitles(hPtAsymCorr[c][p],"Corrected Asymmetry","","events",false);
      if( c == 0 ) hPtAsymCorr[c][p]->SetMarkerStyle(20);
      else hPtAsymCorr[c][p]->SetFillColor(mcColor);
      hPtAsymCorrWidthStd[c]->SetBinContent(1+p,hPtAsymCorr[c][p]->GetRMS());
      hPtAsymCorrWidthStd[c]->SetBinError(1+p,hPtAsymCorr[c][p]->GetRMSError());
      if( hPtAsymCorr[c][p]->Fit("gaus","0QIR","",hPtAsymCorr[c][p]->GetMean()-2.*hPtAsymCorr[c][p]->GetRMS(),hPtAsymCorr[c][p]->GetMean()+2.*hPtAsymCorr[c][p]->GetRMS()) == 0 ) {
	hPtAsymCorrWidthGauss[c]->SetBinContent(1+p,std::abs(hPtAsymCorr[c][p]->GetFunction("gaus")->GetParameter(2)));
	hPtAsymCorrWidthGauss[c]->SetBinError(1+p,hPtAsymCorr[c][p]->GetFunction("gaus")->GetParError(2));
      }
    }
  }
  TH1 *hPtAsymCorrWidthRatioGauss = util::HistOps::createRatioPlot(hPtAsymCorrWidthGauss[0],hPtAsymCorrWidthGauss[1],"#sigma(Asymmetry) Data / MC",0.7,1.5);
  TH1 *hPtAsymCorrWidthRatioStd = util::HistOps::createRatioPlot(hPtAsymCorrWidthStd[0],hPtAsymCorrWidthStd[1],"StdDev(Asymmetry) Data / MC",0.7,1.5);



  // Labels
  std::vector<TPaveText*> labPtAveBin(ptAveBinEdges.size()-1);
  for(size_t i = 0; i < labPtAveBin.size(); ++i) {
    labPtAveBin[i] = util::LabelFactory::createPaveText(3,-0.7);
    labPtAveBin[i]->AddText("L = "+util::toTString(lumi)+" pb^{-1},  "+util::toTString(minEta)+" < |#eta| < "+util::toTString(maxEta));
    labPtAveBin[i]->AddText("|#Delta#phi| > "+util::toTString(minDeltaPhi)+",  p^{rel}_{T} < "+util::toTString(maxPt3Rel));
    labPtAveBin[i]->AddText(util::toTString(ptAveBinEdges[i])+" < p^{ave}_{T} < "+util::toTString(ptAveBinEdges[i+1])+" GeV");
  }
  TPaveText *labSpec = util::LabelFactory::createPaveText(3,-0.7);
  labSpec->AddText("L = "+util::toTString(lumi)+" pb^{-1},  "+util::toTString(minEta)+" < |#eta| < "+util::toTString(maxEta));
  labSpec->AddText("|#Delta#phi| > "+util::toTString(minDeltaPhi)+",  p^{rel}_{T} < "+util::toTString(maxPt3Rel));
  labSpec->AddText("p^{ave}_{T} > "+util::toTString(minPtAve)+" GeV");

  TLegend *leg = util::LabelFactory::createLegendCol(2,0.3);
  leg->AddEntry(hPtCorr[0][0],"Data","P");
  leg->AddEntry(hPtCorr[1][0],"MC","F");



  // Plots
  TCanvas *can1 = new TCanvas("canPtAveCorr","PtAve Corr",500,500);
  can1->cd();
  util::HistOps::setYRange(hPtAveCorr[1],3);
  hPtAveCorr[1]->Draw("HISTE");
  hPtAveCorr[0]->Draw("PE1same");
  labSpec->Draw("same");
  leg->Draw("same");
  can1->SaveAs(prefix+"PtAveCorr.eps","eps");
  hPtAveCorr[1]->GetYaxis()->SetRangeUser(3E-2,7*pow(10,log10(hPtAveCorr[1]->GetMaximum())+3));
  hPtAveCorr[1]->Draw("HISTE");
  hPtAveCorr[0]->Draw("PE1same");
  can1->SetLogy();
  can1->SaveAs(prefix+"PtAveCorrLogy.eps","eps");

  TCanvas *can2 = new TCanvas("canPtAveCorrLog","PtAve Corr (Log)",500,500);
  can2->cd();
  util::HistOps::setYRange(hPtAveCorrLog[1],3);
  hPtAveCorrLog[1]->Draw("HISTE");
  hPtAveCorrLog[0]->Draw("PE1same");
  labSpec->Draw("same");
  leg->Draw("same");
  can2->SetLogx();
  can2->SaveAs(prefix+"PtAveCorrLogx.eps","eps");
  hPtAveCorrLog[1]->GetYaxis()->SetRangeUser(3E-2,7*pow(10,log10(hPtAveCorrLog[1]->GetMaximum())+3));
  hPtAveCorrLog[1]->Draw("HISTE");
  hPtAveCorrLog[0]->Draw("PE1same");
  labSpec->Draw("same");
  leg->Draw("same");
  can2->SetLogy();
  can2->SaveAs(prefix+"PtAveCorrLogxy.eps","eps");

  TCanvas *can3 = new TCanvas("canDeltaPhi","Delta Phi",500,500);
  can3->cd();
  hDeltaPhi[1]->Draw("HISTE");
  hDeltaPhi[0]->Draw("PE1same");
  can3->SaveAs(prefix+"DeltaPhi.eps","eps");

  TCanvas *can4 = new TCanvas("canEta","Eta",500,500);
  can4->cd();
  hEta[1]->Draw("HISTE");
  hEta[0]->Draw("PE1same");
  can4->SaveAs(prefix+"Eta.eps","eps");

  TCanvas *can5 = new TCanvas("canPtAsymCorrWidthGauss","PtAsym Width Gauss",500,500);
  can5->cd();
  util::HistOps::setYRange(hPtAsymCorrWidthGauss[1],3);
  hPtAsymCorrWidthGauss[1]->Draw("HISTE");
  hPtAsymCorrWidthGauss[0]->Draw("PE1same");
  labSpec->Draw("same");
  leg->Draw("same");
  can5->SetLogx();
  can5->SaveAs(prefix+"PtAsymCorrWidthGauss.eps","eps");

  TCanvas *can6 = new TCanvas("canPtAsymCorrWidthGaussRatio","PtAsym Width Gauss Ratio",500,500);
  can6->cd();
  util::HistOps::createRatioFrame(hPtAsymCorrWidthRatioGauss,"#sigma(Asymmetry) Data / MC",0.7,1.8)->Draw();
  hPtAsymCorrWidthRatioGauss->Draw("PE1same");
  labSpec->Draw("same");
  leg->Draw("same");
  can6->SetLogx();
  can6->SaveAs(prefix+"PtAsymCorrWidthGaussRatio.eps","eps");

  TCanvas *can7 = new TCanvas("canPtAsymCorrWidthStd","PtAsym Width Std",500,500);
  can7->cd();
  util::HistOps::setYRange(hPtAsymCorrWidthStd[1],3);
  hPtAsymCorrWidthStd[1]->Draw("HISTE");
  hPtAsymCorrWidthStd[0]->Draw("PE1same");
  labSpec->Draw("same");
  leg->Draw("same");
  can7->SetLogx();
  can7->SaveAs(prefix+"PtAsymCorrWidthStd.eps","eps");

  TCanvas *can8 = new TCanvas("canPtAsymCorrWidthStdRatio","PtAsym Width Std Ratio",500,500);
  can8->cd();
  util::HistOps::createRatioFrame(hPtAsymCorrWidthRatioStd,"StdDev(Asymmetry) Data / MC",0.7,1.8)->Draw();
  hPtAsymCorrWidthRatioStd->Draw("PE1same");
  labSpec->Draw("same");
  leg->Draw("same");
  can8->SetLogx();
  can8->SaveAs(prefix+"PtAsymCorrWidthStdRatio.eps","eps");

  TCanvas *can9 = new TCanvas("canPtJet1vs2Uncorr","1 vs 2 uncorr",500,500);
  can9->cd();
  util::HistOps::setMarginsColz(can9);
  hPtJet1vs2Uncorr[0]->Draw("COLZ");
  can9->SaveAs(prefix+"PtJet1vs2Uncorr.eps","eps");
  hPtJet1vs2UncorrLog[0]->Draw("COLZ");
  can9->SetLogx();
  can9->SetLogy();
  can9->SaveAs(prefix+"PtJet1vs2UncorrLog.eps","eps");

  TCanvas *can10 = new TCanvas("canPtJet1vs2Corr","1 vs 2 corr",500,500);
  can10->cd();
  util::HistOps::setMarginsColz(can10);
  hPtJet1vs2Corr[0]->Draw("COLZ");
  can10->SaveAs(prefix+"PtJet1vs2Corr.eps","eps");
  hPtJet1vs2CorrLog[0]->Draw("COLZ");
  can10->SetLogx();
  can10->SetLogy();
  can10->SaveAs(prefix+"PtJet1vs2CorrLog.eps","eps");

  TCanvas *can11 = new TCanvas("canPtAsymVsEMF","Asym vs emf",500,500);
  can11->cd();
  util::HistOps::setMarginsColz(can11);
  hPtAsymVsEMF[0]->Draw("COLZ");
  can11->SetLogz();
  can11->SaveAs(prefix+"PtAsymUncorrVsEMF.eps","eps");

  for(size_t i = 0; i < hPtUncorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtUncorr"+util::toTString(i),"Pt Uncorr "+util::toTString(1+i),500,500);
    can->cd();
    hPtUncorr[1][i]->GetYaxis()->SetRangeUser(3E-2,7*pow(10,log10(hPtUncorr[1][i]->GetMaximum())+3));
    hPtUncorr[1][i]->Draw("HISTE");
    hPtUncorr[0][i]->Draw("PE1same");
    labSpec->Draw("same");
    leg->Draw("same");
    can->SetLogy();
    can->SaveAs(prefix+"PtUncorr_Jet"+util::toTString(1+i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtCorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtCorr"+util::toTString(i),"Pt Corr "+util::toTString(1+i),500,500);
    can->cd();
    hPtCorr[1][i]->GetYaxis()->SetRangeUser(3E-2,7*pow(10,log10(hPtCorr[1][i]->GetMaximum())+3));
    hPtCorr[1][i]->Draw("HISTE");
    hPtCorr[0][i]->Draw("PE1same");
    labSpec->Draw("same");
    leg->Draw("same");
    can->SetLogy();
    can->SaveAs(prefix+"PtCorr_Jet"+util::toTString(1+i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtAsymUncorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtAsym"+util::toTString(i),"Asym Uncorr "+util::toTString(1+i),500,500);
    can->cd();
    util::HistOps::setYRange(hPtAsymUncorr[1][i],3);
    hPtAsymUncorr[1][i]->Draw("HISTE");
    hPtAsymUncorr[0][i]->Draw("PE1same");
    labPtAveBin[i]->Draw("same");
    leg->Draw("same");
    can->SaveAs(prefix+"PtAsymUncorr_PtAveBin"+util::toTString(i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtAsymUncorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtAsymLog"+util::toTString(i),"Asym Uncorr "+util::toTString(1+i)+" (Log)",500,500);
    can->cd();
    hPtAsymUncorr[1][i]->GetYaxis()->SetRangeUser(3E-2,7*pow(10,log10(hPtAsymUncorr[1][i]->GetMaximum())+3));
    hPtAsymUncorr[1][i]->Draw("HISTE");
    hPtAsymUncorr[0][i]->Draw("PE1same");
    labPtAveBin[i]->Draw("same");
    leg->Draw("same");
    can->SetLogy();
    can->SaveAs(prefix+"PtAsymUncorrLog_PtAveBin"+util::toTString(i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtAsymCorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtAsymCorr"+util::toTString(i),"Asym Corr "+util::toTString(1+i),500,500);
    can->cd();
    util::HistOps::setYRange(hPtAsymCorr[1][i],3);
    hPtAsymCorr[1][i]->Draw("HISTE");
    hPtAsymCorr[0][i]->Draw("PE1same");
    labPtAveBin[i]->Draw("same");
    leg->Draw("same");
    can->SaveAs(prefix+"PtAsymCorr_PtAveBin"+util::toTString(i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtAsymCorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtAsymCorrLog"+util::toTString(i),"Asym Corr "+util::toTString(1+i)+" (Log)",500,500);
    can->cd();
    hPtAsymCorr[1][i]->GetYaxis()->SetRangeUser(3E-2,7*pow(10,log10(hPtAsymCorr[1][i]->GetMaximum())+3));
    hPtAsymCorr[1][i]->Draw("HISTE");
    hPtAsymCorr[0][i]->Draw("PE1same");
    labPtAveBin[i]->Draw("same");
    leg->Draw("same");
    can->SetLogy();
    can->SaveAs(prefix+"PtAsymCorrLog_PtAveBin"+util::toTString(i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtBiasAsymCorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtBiasAsymCorr"+util::toTString(i),"Biased Asym Corr "+util::toTString(1+i),500,500);
    can->cd();
    util::HistOps::setYRange(hPtBiasAsymCorr[1][i],3);
    hPtBiasAsymCorr[1][i]->Draw("HISTE");
    hPtBiasAsymCorr[0][i]->Draw("PE1same");
    labPtAveBin[i]->Draw("same");
    leg->Draw("same");
    can->SaveAs(prefix+"PtBiasAsymCorr_PtAveBin"+util::toTString(i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtBiasAsymCorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtBiasAsymCorrLog"+util::toTString(i),"Biased Asym Corr "+util::toTString(1+i)+" (Log)",500,500);
    can->cd();
    hPtBiasAsymCorr[1][i]->GetYaxis()->SetRangeUser(3E-2,7*pow(10,log10(hPtBiasAsymCorr[1][i]->GetMaximum())+3));
    hPtBiasAsymCorr[1][i]->Draw("HISTE");
    hPtBiasAsymCorr[0][i]->Draw("PE1same");
    labPtAveBin[i]->Draw("same");
    leg->Draw("same");
    can->SetLogy();
    can->SaveAs(prefix+"PtBiasAsymCorrLog_PtAveBin"+util::toTString(i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtAsymVsPtJet3Uncorr.size(); ++i) {
    TCanvas *can = new TCanvas("canPtAsymVsPtJet3Uncorr"+util::toTString(i),"Asym vs pt3",500,500);
    can->cd();
    hPtAsymVsPtJet3Uncorr[i]->Draw("COLZ");
    can->SetLogx();
    can->SaveAs(prefix+"PtAsymVsPtJet3_"+(i==0 ? "Data" : "MC" )+".eps","eps");
  }


  // Cut-Flow
  std::cout << std::endl;
  std::cout << "Total   : " << nTotal[0] << std::endl;
  std::cout << "PtAve   : " << (nTotal[0] -= nPtAve[0]) << std::endl;
  std::cout << "Eta     : " << (nTotal[0] -= nEta[0]) << std::endl;
  std::cout << "DPhi    : " << (nTotal[0] -= nDeltaPhi[0]) << std::endl;
  std::cout << "Pt3     : " << (nTotal[0] -= nPt3Rel[0]) << std::endl;
  std::cout << "JetID   : " << (nTotal[0] -= nJetID[0]) << std::endl;
}

