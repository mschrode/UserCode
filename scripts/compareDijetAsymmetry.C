#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include "TString.h"
#include "TVector2.h"

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
  util::StyleSettings::presentation();

  const TString prefix = "JetMET_Run2010A-PromptReco-v4_614nb_Pt3RelCuts_";
  const int maxNJet = 50;
  const double minEMF = 0.01;
  const int minJetN90Hits = 2;
  const double maxJetFHPD = 0.98;
  const double maxJetFRBX = 0.98;
  const double minDeltaPhi = 2.7;
  const double maxEta = 1.;
  const double minPtAve = 20.;
  const double maxPt3Rel = 0.1;
  const double lumi = 0.614;
  const double lumiFrac = lumi/0.1;
  TRandom *rand = new TRandom(0);

  // 0: Data
  // 1: MC
  unsigned int nChains = 2;

  //TChain *data = createTChain("/afs/naf.desy.de/user/m/mschrode/Kalibri/input/Kalibri_JetMETTau_Run2010A-Jul15thReRec-v1A"); // 70nb
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
  std::vector<TH1*> hPtAveCorr(nChains);
  std::vector<TH1*> hPtAveCorrLog(nChains);
  std::vector<TH2*> hPtAsymVsPtAveUncorr(nChains);
  std::vector<TH2*> hPtAsymVsPtAveCorr(nChains);
  std::vector<TH2*> hPtAsymVsPtJet3Uncorr(nChains);
  std::vector<TH1*> hDeltaPhi(nChains);
  std::vector<TH1*> hEta(nChains);
  for(size_t c = 0; c < nChains; ++c) {
    TH1 *h = 0;
    for(int j = 0; j < 3; ++j) {
      h = util::HistOps::createTH1D("hPtUncorr"+util::toTString(c)+"_"+util::toTString(j+1),50,ptAveBinEdges.front(),ptAveBinEdges.back(),"p_{T,"+util::toTString(j+1)+"}","GeV","jets",false);
      h->Sumw2();
      if( c == 0 ) h->SetMarkerStyle(20);
      hPtUncorr[c].push_back(h);

      h = util::HistOps::createTH1D("hPtCorr"+util::toTString(c)+"_"+util::toTString(j+1),50,ptAveBinEdges.front(),ptAveBinEdges.back(),"Corrected p_{T,"+util::toTString(j+1)+"}","GeV","jets",false);
      h->Sumw2();
      if( c == 0 ) h->SetMarkerStyle(20);
      hPtCorr[c].push_back(h);
    }
    h = util::HistOps::createTH1D("hPtAveCorr"+util::toTString(c),50,ptAveBinEdges.front(),ptAveBinEdges.back(),"Corrected p^{ave}_{T}","GeV","jets",false);
    h->Sumw2();
    if( c == 0 ) h->SetMarkerStyle(20);
    hPtAveCorr[c] = h;

    h = util::HistOps::createTH1D("hPtAveCorrLog"+util::toTString(c),ptAveBinEdgesLog.size()-1,&(ptAveBinEdgesLog.front()),"Corrected p^{ave}_{T}","GeV","jets",false);
    h->Sumw2();
    if( c == 0 ) h->SetMarkerStyle(20);
    hPtAveCorrLog[c] = h;

    h = util::HistOps::createTH1D("hDeltaPhi"+util::toTString(c),50,0.,M_PI,"|#Delta#phi_{12}|","","events",false);
    h->Sumw2();
    if( c == 0 ) h->SetMarkerStyle(20);
    hDeltaPhi[c] = h;

    h = util::HistOps::createTH1D("hEta"+util::toTString(c),51,-5.1,5.1,"#eta_{1,2}","","jets",false);
    h->Sumw2();
    if( c == 0 ) h->SetMarkerStyle(20);
    hEta[c] = h;

    TH2 *h2 = new TH2D("hPtAsymVsPtAveUncorr"+util::toTString(c),"",ptAveBinEdges.size()-1,&(ptAveBinEdges.front()),31,-1.,1.);
    h2->Sumw2();
    hPtAsymVsPtAveUncorr[c] = h2;

    h2 = new TH2D("hPtAsymVsPtAveCorr"+util::toTString(c),"",ptAveBinEdges.size()-1,&(ptAveBinEdges.front()),31,-1.,1.);
    h2->Sumw2();
    hPtAsymVsPtAveCorr[c] = h2;

    h2 = new TH2D("hPtAsymVsPtJet3Uncorr"+util::toTString(c),"p_{T,3} (GeV);Asymmetry",ptJet3BinEdgesLog.size()-1,&(ptJet3BinEdgesLog.front()),31,-1.,1.);
    h2->Sumw2();
    hPtAsymVsPtJet3Uncorr[c] = h2;
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
	else if( std::abs(jetEta[0]) > maxEta || std::abs(jetEta[1]) > maxEta ) isGoodEvt = false;
	else if( std::abs(TVector2::Phi_mpi_pi(jetPhi[0]-jetPhi[1])) < minDeltaPhi ) isGoodEvt = false;
	else if( (nObjJet>2 ? jetPt[2] : 0.) > maxPt3Rel*ptAveCorr ) isGoodEvt = false;
	else if( jetEMF[0] < minEMF || jetEMF[1] < minEMF ) isGoodEvt = false;
	else if( jetN90Hits[0] < minJetN90Hits || jetN90Hits[1] < minJetN90Hits ) isGoodEvt = false;
	else if( jetFHPD[0] > maxJetFHPD || jetFHPD[1] > maxJetFHPD ) isGoodEvt = false;
	else if( jetFRBX[0] > maxJetFRBX || jetFRBX[1] > maxJetFRBX ) isGoodEvt = false;

	if( isGoodEvt ) {
	  hPtUncorr[c][0]->Fill(jetPt[0],weight);
	  hPtUncorr[c][1]->Fill(jetPt[1],weight);
	  if( nObjJet > 2 ) hPtUncorr[c][2]->Fill(jetPt[2],weight);

	  double ptAsym = jetPt[0]+jetPt[1];
	  if( ptAsym > 0. ) {
	    ptAsym = (jetPt[0]-jetPt[1])/ptAsym;
	    if( rand->Uniform()>0.5 ) ptAsym *= -1.;
	    hPtAsymVsPtAveUncorr[c]->Fill(0.5*(jetPt[0]+jetPt[1]),ptAsym,weight);
	    if( nObjJet > 2 ) hPtAsymVsPtJet3Uncorr[c]->Fill(jetPt[2],ptAsym,weight);
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
  	else if( std::abs(jetEta[jIdx[0]->idx_]) > maxEta || std::abs(jetEta[jIdx[1]->idx_]) > maxEta ) {
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
	else if( jetEMF[jIdx[0]->idx_] < minEMF || jetEMF[jIdx[1]->idx_] < minEMF ) {
	  nJetID[c]++;
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
	else if( jetFRBX[jIdx[0]->idx_] > maxJetFRBX || jetFRBX[jIdx[1]->idx_] > maxJetFRBX ) {
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
	  hEta[c]->Fill(std::abs(jetEta[jIdx[0]->idx_]),weight);
	  hEta[c]->Fill(std::abs(jetEta[jIdx[1]->idx_]),weight);

  	  double ptAsym = jIdx[0]->pt_+jIdx[1]->pt_;
  	  if( ptAsym > 0. ) {
  	    ptAsym = (jIdx[0]->pt_-jIdx[1]->pt_)/ptAsym;
  	    if( rand->Uniform()>0.5 ) ptAsym *= -1.;
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
  for(unsigned int c = 0; c < nChains; ++c) {
    util::HistOps::fillSlices(hPtAsymVsPtAveUncorr[c],hPtAsymUncorr[c],"hPtAsymUncorr"+util::toTString(c));
    util::HistOps::fillSlices(hPtAsymVsPtAveCorr[c],hPtAsymCorr[c],"hPtAsymCorr"+util::toTString(c));
    for(size_t p = 0; p < hPtAsymUncorr[c].size(); ++p) {
      util::HistOps::setAxisTitles(hPtAsymUncorr[c][p],"Asymmetry","","events",false);
      if( c == 0 ) hPtAsymUncorr[c][p]->SetMarkerStyle(20);
      hPtAsymUncorr[c][p]->SetTitle(util::toTString(ptAveBinEdges[p])+" < p^{ave}_{T} < "+util::toTString(ptAveBinEdges[p+1])+" GeV");

      util::HistOps::setAxisTitles(hPtAsymCorr[c][p],"Corrected Asymmetry","","events",false);
      if( c == 0 ) hPtAsymCorr[c][p]->SetMarkerStyle(20);
      hPtAsymCorr[c][p]->SetTitle(util::toTString(ptAveBinEdges[p])+" < p^{ave}_{T} < "+util::toTString(ptAveBinEdges[p+1])+" GeV");
    }
  }


  // Plots
  TCanvas *can1 = new TCanvas("canPtAveCorr","PtAve Corr",500,500);
  can1->cd();
  hPtAveCorr[1]->Draw("HISTE");
  hPtAveCorr[0]->Draw("PE1same");
  can1->SaveAs(prefix+"PtAveCorr.eps","eps");
  can1->SetLogy();
  can1->SaveAs(prefix+"PtAveCorrLogy.eps","eps");

  TCanvas *can2 = new TCanvas("canPtAveCorrLog","PtAve Corr (Log)",500,500);
  can2->cd();
  hPtAveCorrLog[1]->Draw("HISTE");
  hPtAveCorrLog[0]->Draw("PE1same");
  can2->SetLogx();
  can2->SaveAs(prefix+"PtAveCorrLogx.eps","eps");
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

  for(size_t i = 0; i < hPtUncorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtUncorr"+util::toTString(i),"Pt Uncorr "+util::toTString(1+i),500,500);
    can->cd();
    hPtUncorr[1][i]->Draw("HISTE");
    hPtUncorr[0][i]->Draw("PE1same");
    can->SetLogy();
    can->SaveAs(prefix+"PtUncorr_Jet"+util::toTString(1+i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtCorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtCorr"+util::toTString(i),"Pt Corr "+util::toTString(1+i),500,500);
    can->cd();
    hPtCorr[1][i]->Draw("HISTE");
    hPtCorr[0][i]->Draw("PE1same");
    can->SetLogy();
    can->SaveAs(prefix+"PtCorr_Jet"+util::toTString(1+i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtAsymUncorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtAsym"+util::toTString(i),"Asym Uncorr "+util::toTString(1+i),500,500);
    can->cd();
    hPtAsymUncorr[1][i]->Draw("HISTE");
    hPtAsymUncorr[0][i]->Draw("PE1same");
    can->SaveAs(prefix+"PtAsymUncorr_PtAveBin"+util::toTString(i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtAsymUncorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtAsymLog"+util::toTString(i),"Asym Uncorr "+util::toTString(1+i)+" (Log)",500,500);
    can->cd();
    hPtAsymUncorr[1][i]->Draw("HISTE");
    hPtAsymUncorr[0][i]->Draw("PE1same");
    can->SetLogy();
    can->SaveAs(prefix+"PtAsymUncorrLog_PtAveBin"+util::toTString(i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtAsymCorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtAsymCorr"+util::toTString(i),"Asym Corr "+util::toTString(1+i),500,500);
    can->cd();
    hPtAsymCorr[1][i]->Draw("HISTE");
    hPtAsymCorr[0][i]->Draw("PE1same");
    can->SaveAs(prefix+"PtAsymCorr_PtAveBin"+util::toTString(i)+".eps","eps");
  }
  for(size_t i = 0; i < hPtAsymCorr[0].size(); ++i) {
    TCanvas *can = new TCanvas("canPtAsymCorrLog"+util::toTString(i),"Asym Corr "+util::toTString(1+i)+" (Log)",500,500);
    can->cd();
    hPtAsymCorr[1][i]->Draw("HISTE");
    hPtAsymCorr[0][i]->Draw("PE1same");
    can->SetLogy();
    can->SaveAs(prefix+"PtAsymCorrLog_PtAveBin"+util::toTString(i)+".eps","eps");
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

