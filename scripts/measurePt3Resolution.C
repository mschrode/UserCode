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
#include "TString.h"
#include "TVector2.h"

#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/utils.h"
#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/HistOps.h"
#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/LabelFactory.h"
#include "/afs/naf.desy.de/user/m/mschrode/UserCode/mschrode/util/StyleSettings.h"



TString prefix_ = "Test3rdJetRes";
const int maxNJet_ = 50;
const double minDeltaPhi_ = 150.;
const double maxDeltaPhi_ = 170.;
const double minEta_ = 0.;
const double maxEta_ = 1.3;
const double minPtAve_ = 60.;
const double maxPt4_ = 10.;
const double ptRecoilMin_ = 20.;
const double ptRecoilMax_ = 60.;
const double lumi_ = 2.9;
const double lumiMC_ = 0.1;
const double resPar_[3] = { 3.04407, 1.16007, 0.0348195 }; // MC truth Eta 0. - 1.3


TFile *outFile_ = 0;
double minDeltaPhiRad_;
double maxDeltaPhiRad_;
double lumiFrac_;
int mcColor_;
TF1 *res_ = 0;
TH1 *hPt3Gen_ = 0;
TH1 *hPtRecoil_ = 0;
TH1 *hPtAveGen_ = 0;
TH1 *hPtAve_ = 0;
TH2 *h2Pt3OverPtRecoilVsPtRecoil_ = 0;
TH2 *h2Pt3OverPtRecoilVsPt3Gen_ = 0;
TH2 *h2MCResVsPtRecoil_ = 0;
TH2 *h2MCResVsPt3Gen_ = 0;


void init() {
  std::cout << "Init... " << std::flush;

  minDeltaPhiRad_ = M_PI*minDeltaPhi_/180.;
  maxDeltaPhiRad_ = M_PI*maxDeltaPhi_/180.;
  lumiFrac_ = lumi_/lumiMC_;
  mcColor_ = 5;
  res_ = new TF1("res","sqrt([0]*[0]/x/x + [1]*[1]/x + [2]*[2])",1.,1000.);
  for(int i = 0; i < 3; ++i) {
    res_->SetParameter(i,resPar_[i]);
  }

  outFile_ = new TFile(prefix_+".root","RECREATE");
  prefix_ += "_";

  hPt3Gen_ = util::HistOps::createTH1D("hPt3Gen_",30,0.,ptRecoilMax_,"p_{T}","GeV","events");

  hPtRecoil_ = static_cast<TH1D*>(hPt3Gen_->Clone("hPtRecoil_"));
  hPtRecoil_->SetMarkerStyle(20);

  hPtAveGen_ = util::HistOps::createTH1D("hPtAveGen_",30,10.,200.,"p^{ave}_{T}","GeV","events");

  hPtAve_ = static_cast<TH1D*>(hPtAveGen_->Clone("hPtAve_"));
  hPtAve_->SetMarkerStyle(20);

  h2Pt3OverPtRecoilVsPtRecoil_ = new TH2D("h2Pt3OverPtRecoilVsPtRecoil_",
					  ";p^{recoil}_{T} (GeV);p_{T,3} / p^{recoil}_{T}",
					  2,ptRecoilMin_,ptRecoilMax_,21,0.,2.);

  h2Pt3OverPtRecoilVsPt3Gen_ = static_cast<TH2D*>(h2Pt3OverPtRecoilVsPtRecoil_->Clone("h2Pt3OverPtRecoilVsPt3Gen_"));
  h2Pt3OverPtRecoilVsPt3Gen_->SetTitle(";p^{gen}_{T,3} (GeV);p_{T,3} / p^{recoil}_{T}");

  h2MCResVsPtRecoil_ = static_cast<TH2D*>(h2Pt3OverPtRecoilVsPtRecoil_->Clone("h2MCResVsPtRecoil_"));
  h2MCResVsPtRecoil_->SetTitle(";p^{recoil}_{T,3} (GeV);p_{T,3} / p^{gen}_{T,3}");

  h2MCResVsPt3Gen_ = static_cast<TH2D*>(h2Pt3OverPtRecoilVsPtRecoil_->Clone("h2MCResVsPt3Gen_"));
  h2MCResVsPt3Gen_->SetTitle(";p^{gen}_{T,3} (GeV);p_{T,3} / p^{gen}_{T,3}");


  std::cout << "ok" << std::endl;
}


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



void fillHistograms(int nMaxEvts = -1, int prescale = 1) {

  TChain *chain = createTChain("/afs/naf.desy.de/user/m/mschrode/Kalibri/input/Kalibri_Summer10QCDDiJet_AK5Calo_PtHat30");

  if( chain ) {
    // Read variables
    std::vector<JetIndex*> jIdx;
    jIdx.resize(8,0);
    // Init read quantities
    float weight = 1.;
    int nObjJet = 0;
    float jetPt[maxNJet_];
    float genJetPt[maxNJet_];
    float jetEta[maxNJet_];
    float jetPhi[maxNJet_];
    float jetCorrL2L3[maxNJet_];
    bool jetID[maxNJet_];

    // Set branch addresses
    chain->SetBranchAddress("Weight",&weight);
    chain->SetBranchAddress("NobjJet",&nObjJet);
    chain->SetBranchAddress("JetPt",jetPt);
    chain->SetBranchAddress("GenJetPt",genJetPt);
    chain->SetBranchAddress("JetEta",jetEta);
    chain->SetBranchAddress("JetPhi",jetPhi);
    chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
    chain->SetBranchAddress("JetIDLoose",jetID);
      
    // Loop over tree entries and fill histograms
    int nEntries = chain->GetEntries();
    if( nMaxEvts > 0 && nEntries > nMaxEvts ) nEntries = nMaxEvts;
    for(int n = 0; n < nEntries; n = n + prescale) {
      if( n%10000 == 0 ) std::cout << " Entry " << n << std::endl;
      chain->GetEntry(n);

      if( nObjJet > maxNJet_ ) {
	std::cerr << "WARNING: nObjJet = " << nObjJet << " > " << maxNJet_ << ". Skipping event.\n";
	continue;
      }
      
      if( nObjJet < 3 ) continue;

      bool isGoodEvt = true;
      
      // Sort by corrected pt
      unsigned int nJets = nObjJet;
      if( jIdx.size() < nJets ) nJets = jIdx.size();
      for(size_t i = 0; i < nJets; ++i) {
	jIdx[i] = new JetIndex(i,jetPt[i]*jetCorrL2L3[i]);
      }
      std::sort(jIdx.begin(),jIdx.begin()+nJets,JetIndex::ptGreaterThan);
      for(size_t i = 0; i < nJets-1; ++i) {
	//	std::cout << "Pt(" << i << ") = " << jIdx[i]->pt_ << std::endl;
 	if( jIdx[i+1]->pt_ > jIdx[i]->pt_ )
 	  std::cerr << "ERROR: Pt(" << i+1 << ") > Pt(" << i << ")" << std::endl;
      }

      double ptAve = 0.5*(jIdx[0]->pt_+jIdx[1]->pt_);
      double deltaPhi = std::abs(TVector2::Phi_mpi_pi(jetPhi[jIdx[0]->idx_]-jetPhi[jIdx[1]->idx_]));
      double pt4 = ( nJets > 3 ? jIdx[3]->pt_ : 0. );

      // Select events
      if( ptAve < minPtAve_ ) isGoodEvt = false;
      else if( std::abs(jetEta[jIdx[0]->idx_]) < minEta_ || std::abs(jetEta[jIdx[0]->idx_]) > maxEta_ ||
	       std::abs(jetEta[jIdx[1]->idx_]) < minEta_ || std::abs(jetEta[jIdx[1]->idx_]) > maxEta_ ||
	       std::abs(jetEta[jIdx[2]->idx_]) < minEta_ || std::abs(jetEta[jIdx[2]->idx_]) > maxEta_)
	isGoodEvt = false;
      else if( deltaPhi < minDeltaPhiRad_ || deltaPhi > maxDeltaPhiRad_) isGoodEvt = false;
      else if( pt4 > maxPt4_ ) isGoodEvt = false;
      else if( !(jetID[jIdx[0]->idx_] && jetID[jIdx[1]->idx_] && jetID[jIdx[2]->idx_]) ) isGoodEvt = false;
      
      // Fill histograms
      if( isGoodEvt ) {
	double ptRecoil = std::abs(cos(deltaPhi/2.))*(jIdx[0]->pt_+jIdx[1]->pt_);
	double pt3 = jIdx[2]->pt_;
	double pt3Gen = genJetPt[jIdx[2]->idx_];

	hPt3Gen_->Fill(pt3Gen,weight);
	hPtRecoil_->Fill(ptRecoil,weight);
	hPtAve_->Fill(ptAve,weight);
	hPtAveGen_->Fill(0.5*(genJetPt[jIdx[0]->idx_]+genJetPt[jIdx[1]->idx_]),weight);

	h2Pt3OverPtRecoilVsPtRecoil_->Fill( ptRecoil, pt3/ptRecoil, weight );
	h2Pt3OverPtRecoilVsPt3Gen_->Fill(   pt3Gen,   pt3/ptRecoil, weight );
	h2MCResVsPtRecoil_->Fill(           ptRecoil, pt3/pt3Gen,   weight );
	h2MCResVsPt3Gen_->Fill(             pt3Gen,   pt3/pt3Gen,   weight );
      }
      for(size_t i = 0; i < nJets; ++i) {
	delete jIdx[i];
      }
    } // End of loop over entries
  }
}


void fitRecoilResolution(const TH2 *h2, util::HistVec &hResponse, TH1 *hResolution, std::vector<TF1*> &fits) {
  if( hResolution ) {
    if( hResolution->GetNbinsX() == static_cast<int>(hResponse.size()) ) {
      fits.clear();
      for(size_t i = 0; i < hResponse.size(); ++i) {
	double aveDeltaPhi = 0.5*(minDeltaPhiRad_+maxDeltaPhiRad_);
	// Average PtRecoil and ptScale
	double ptAveRecoil = h2->GetXaxis()->GetBinCenter(1+i);
	double ptScale = ptAveRecoil / 2. / std::abs(cos(aveDeltaPhi/2.));
	// PtRecoil resolution / ptAveRecoil from higher pt resolution
	double sPtRecoil = res_->Eval(ptScale)/sqrt(2.);

	// Fit pt3/ptRecoil
	TString name = hResponse[i]->GetName();
	name += "_fit";
	fits.push_back(new TF1(name,"gaus",0.,2.));
	hResponse[i]->Fit(fits[i],"0QI");
	double sPt3OverPtRecoil = fits[i]->GetParameter(2);
	double tmp = sPt3OverPtRecoil*sPt3OverPtRecoil - sPtRecoil*sPtRecoil;
	double sPt3 = ( tmp > 0. ? sqrt(tmp) : 1. );
	double sPt3Err = fits[i]->GetParError(2)/sPt3;

	std::cout << std::endl;
	std::cout << "ptAveRecoil      = " << ptAveRecoil << std::endl;
	std::cout << "ptScale          = " << ptScale << std::endl;
	std::cout << "sPt3OverPtRecoil = " << sPt3OverPtRecoil << std::endl;
	std::cout << "sPt3             = " << sPt3 << std::endl;

  	hResolution->SetBinContent(1+i,sPt3);
  	hResolution->SetBinError(1+i,sPt3Err);
    }
    } else {
      std::cerr << "ERROR(fitRecoilResolution): incorrect number of bins" << std::endl;
    }
  } else {
    std::cerr << "ERROR(fitRecoilResolution): histogram not initialised" << std::endl;
  }
}

void makePlots() {
  // Control plots
  TCanvas *canPtAve = new TCanvas("canPtAve","PtAve",500,500);
  canPtAve->cd();
  hPtAveGen_->Draw("HISTE");
  hPtAve_->Draw("PE1same");

  TCanvas *canPtRecoil = new TCanvas("canPtRecoil","PtRecoil",500,500);
  canPtRecoil->cd();
  hPt3Gen_->Draw("HISTE");
  hPtRecoil_->Draw("PE1same");


  // Calculate resolution
  util::HistVec hMeasRes_PtRecoil; 
  util::HistVec hMeasRes_PtGen;
  util::HistVec hMCRes_PtRecoil;
  util::HistVec hMCRes_PtGen;
  util::HistOps::fillSlices(h2Pt3OverPtRecoilVsPtRecoil_,hMeasRes_PtRecoil,"hMeasRes_PtRecoil");
  util::HistOps::fillSlices(h2Pt3OverPtRecoilVsPt3Gen_,hMeasRes_PtGen,"hMeasRes_PtGen");
  util::HistOps::fillSlices(h2MCResVsPtRecoil_,hMCRes_PtRecoil,"hMCRes_PtRecoil");
  util::HistOps::fillSlices(h2MCResVsPt3Gen_,hMCRes_PtGen,"hMCRes_PtGen");

  std::vector<TF1*> fMeasResVsPtRecoil;
  std::vector<TF1*> fMeasResVsPtGen;
  std::vector<TF1*> fMCResVsPtRecoilcoil;
  std::vector<TF1*> fMCResVsPtGen;

  TH1 *hMeasResVsPtRecoil = util::HistOps::createTH1D("hMeasResVsPtRecoil",hMeasRes_PtRecoil.size(),ptRecoilMin_,ptRecoilMax_,"p_{T}","GeV","Resolution"); 
  hMeasResVsPtRecoil->SetMarkerStyle(20);
  TH1 *hMeasResVsPtGen = static_cast<TH1D*>(hMeasResVsPtRecoil->Clone("hMeasResVsPtGen"));
  TH1 *hMCResVsPtRecoil = static_cast<TH1D*>(hMeasResVsPtRecoil->Clone("hMCResVsPtRecoil"));
  TH1 *hMCResVsPtGen = static_cast<TH1D*>(hMCResVsPtRecoil->Clone("hMCResVsPtGen"));

  // Resolution of third jet
  fitRecoilResolution(h2Pt3OverPtRecoilVsPtRecoil_,hMeasRes_PtRecoil,hMeasResVsPtRecoil,fMeasResVsPtRecoil);
  
  TCanvas *canRes1 = new TCanvas("canRes1","Resolution 1",500,500);
  canRes1->cd();
  hMeasResVsPtRecoil->Draw("PE1");
  res_->Draw("same");


  for(size_t i = 0; i < hMeasRes_PtRecoil.size(); ++i) {
    TCanvas *can = new TCanvas("canVsPtRecoil"+util::toTString(i),"Resolution vs PtMeas ("+util::toTString(i)+")",500,500);
    can->cd();
    hMCRes_PtRecoil[i]->Draw("HISTE");
    hMeasRes_PtRecoil[i]->SetMarkerStyle(20);
    hMeasRes_PtRecoil[i]->Draw("PE1same");

//     TCanvas *can2 = new TCanvas("canVsPtGen"+util::toTString(i),"Resolution vs PtGen ("+util::toTString(i)+")",500,500);
//     can2->cd();
//     hMCRes_PtGen[i]->Draw("HISTE");
//     hMeasRes_PtGen[i]->SetMarkerStyle(20);
//     hMeasRes_PtGen[i]->Draw("PE1same");
  }

  outFile_->Write();
}


void measurePt3Resolution() {
  util::StyleSettings::presentationNoTitle();

  init();
  fillHistograms();
  makePlots();
}
