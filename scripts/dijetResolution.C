// $Id: $

// Note: This script produces many plots; best use
// root -b when running it!

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TError.h"
#include "TH1D.h"
#include "TH2I.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TVector2.h"


// === Type definitions ===
class Jet {
public:
  Jet(double eGen, double ptGen, double e, double pt, double eta, double phi, double emf, double dR)
    : eGen_(eGen), ptGen_(ptGen), e_(e), pt_(pt), eta_(eta), phi_(phi), emf_(emf), dR_(dR) {};

  bool operator<(const Jet& jet) const { return ( ptGen() < jet.ptGen() ); }

  double eGen() const { return eGen_; }
  double ptGen() const { return ptGen_; }
  double e() const { return e_; }
  double pt() const { return pt_; }
  double eta() const { return eta_; }
  double phi() const { return phi_; }
  double emf() const { return emf_; }
  double dR() const { return dR_; }

private:
  const double eGen_;
  const double ptGen_;
  const double e_;
  const double pt_;
  const double eta_;
  const double phi_;
  const double emf_;
  const double dR_;
};
typedef std::vector<Jet*>::const_iterator JetIt;


class Event {
public:
  Event() {
    isSelected_ = false;
    isDijet_ = false;
    ht_ = 0.;
    htGen_ = 0.;
    mht_ = 0.;
    mhtGen_ = 0.;
    mhtDijet_ = 0.;
    mhtNjet_ = 0.;
    maxDeltaR_ = 0.;
  }

  ~Event() {
    for(std::vector<Jet*>::iterator it = jets_.begin(); it != jets_.end(); it++) {
      delete *it;
    }
  }
  
  void addJet(Jet *jet) {
    jets_.push_back(jet);
    init();
  }

  int nJets() const { return static_cast<int>(jets_.size()); }
  JetIt firstJet() const { return jets_.begin(); }
  JetIt lastJet() const { return jets_.end(); }

  bool isDijet() const { return isDijet_; }
  void isDijet(bool isDijet) { isDijet_ = isDijet; }
  bool isSelected() const { return isSelected_; }
  void isSelected(bool isSelected) { isSelected_ = isSelected; }

  double ht() const { return ht_; }
  double htGen() const { return htGen_; }
  double mht() const { return mht_; }
  double mhtGen() const { return mhtGen_; }
  double mhtPredDijet() const { return mhtDijet_; }
  double mhtPredNjet() const { return mhtNjet_; }
  void mhtPredDijet(double mht) { mhtDijet_ = mht; }
  void mhtPredNjet(double mht) { mhtNjet_ = mht; }

  double maxDeltaR() const { return maxDeltaR_; }
  double deltaPhi() const { return nJets() < 2 ? 0.
			      : TVector2::Phi_mpi_pi(jets_[0]->phi() - jets_[1]->phi()); }
  double deltaPhiAbs() const { return std::abs(deltaPhi()); }
  double ptDijet() const { return nJets() < 2 ? 0. : 0.5*(jets_[0]->pt() + jets_[1]->pt()); }
  double rel3rdJetPt() const { return nJets() < 3 ? 0. : jets_[2]->pt() / ptDijet(); }

private:
  std::vector<Jet*> jets_;
  bool isDijet_;
  bool isSelected_;
  double ht_;
  double htGen_;
  double mht_;
  double mhtGen_;
  double mhtDijet_;
  double mhtNjet_;
  double maxDeltaR_;

  void init() {
    std::sort(jets_.begin(),jets_.end());
    ht_ = 0.;
    htGen_ = 0.;
    double mht[2] = { 0., 0. };
    double mhtGen[2] = { 0., 0. };
    maxDeltaR_ = 0.;
    for(JetIt jetIt = firstJet(); jetIt != lastJet(); jetIt++) {
      Jet *jet = *jetIt;

      ht_ += jet->pt();
      htGen_ += jet->ptGen();

      mht[0] += jet->pt()*cos(jet->phi());
      mht[1] += jet->pt()*sin(jet->phi());
      mhtGen[0] += jet->ptGen()*cos(jet->phi());
      mhtGen[1] += jet->ptGen()*sin(jet->phi());

      if( jet->dR() > maxDeltaR_ ) maxDeltaR_ = jet->dR();
    }
    mht_ = sqrt( mht[0]*mht[0] + mht[1]*mht[1] );
    mhtGen_ = sqrt( mhtGen[0]*mhtGen[0] + mhtGen[1]*mhtGen[1] );
  }
};
typedef std::vector<Event*> Data;
typedef std::vector<Event*>::const_iterator DataIt;


class Response {
public:
  Response(const TString &id, const std::vector<double> &eBins, const std::vector<double> &etaBins)
    : eMin_(eBins.front()),
      eMax_(eBins.back()),
      etaMin_(etaBins.front()),
      etaMax_(etaBins.back()) {
    TString name = id;
    name += "::indexContainer";
    indexContainer_ = new TH2I(name,"",
			       eBins.size()-1,&(eBins.front()),
			       etaBins.size()-1,&(etaBins.front()));
    responsePdf_ = std::vector<TH1D*>(nBins());
    for(int etaBin = 0; etaBin < nEtaBins(); etaBin++) {
      for(int eBin = 0; eBin < nEBins(); eBin++) {
	// Global bin index
	int idx = eBin + etaBin*nEBins();
	indexContainer_->SetBinContent(eBin+1,etaBin+1,idx);
	// Create response histogram
	name = id;
	name += "::pdf_";
	name += idx;
	char title[100];
	sprintf(title,"%.0f < E < %.0f GeV, %.1f < #eta < %.1f;Response",
		indexContainer_->GetXaxis()->GetBinLowEdge(eBin+1),
		indexContainer_->GetXaxis()->GetBinUpEdge(eBin+1),
		indexContainer_->GetYaxis()->GetBinLowEdge(etaBin+1),
		indexContainer_->GetYaxis()->GetBinUpEdge(etaBin+1));
	responsePdf_[idx] = new TH1D(name,title,50,0,2);
      }
    }
  }

  ~Response() {
    for(int i = 0; i < nBins(); i++) {
      delete responsePdf_[i];
    }
    delete indexContainer_;
  }

  bool isInRange(double e, double eta) const {
    return e > eMin() && e < eMax() && eta > etaMin() && eta < etaMax();
  }
  double eMin() const { return eMin_; }
  double eMax() const { return eMax_; }
  double etaMin() const { return etaMin_; }
  double etaMax() const { return etaMax_; }
  int nBins() const { return nEBins()*nEtaBins(); }
  int nEBins() const { return indexContainer_->GetNbinsX(); }
  int nEtaBins() const { return indexContainer_->GetNbinsY(); }
  
  double random(double e, double eta) const {
    int idx = bin(e,eta);
    double ran = -1.;
    if( idx >= 0 && idx < nBins() ) {
      if( responsePdf_[idx]->Integral() ) {
	ran = responsePdf_[idx]->GetRandom();
      }
    }
    return ran;
  }
  void fill(double e, double eta, double r) {
    int idx = bin(e,eta);
    if( idx >= 0 && idx < nBins() ) responsePdf_[idx]->Fill(r);
  }
  TH1D *pdf(double e, double eta) {
    int idx = bin(e,eta);
    return idx >= 0 && idx < nBins() ? responsePdf_[idx] : 0;
  }
  TH1D *pdf(int bin) {
    return bin >= 0 && bin < nBins() ? responsePdf_[bin] : 0;
  }

private:
  TH2I *indexContainer_;
  std::vector<TH1D*> responsePdf_;
  double eMin_;
  double eMax_;
  double etaMin_;
  double etaMax_;

  int bin(int eBin, int etaBin) const {
    return static_cast<int>(indexContainer_->GetBinContent(indexContainer_->GetBin(eBin+1,etaBin+1)));
  }
  int bin(double e, double eta) const {
    return static_cast<int>(indexContainer_->GetBinContent(indexContainer_->FindBin(e,eta)));
  }
};


typedef std::vector<TCanvas*> Can;
typedef std::vector<TCanvas*>::iterator CanIt;
typedef std::vector<TH1D*> Hist;
typedef std::vector<TH1D*>::iterator HistIt;



// === Global variables ===
// Auxiliary
const int maxNJet_ = 50;
const int maxNTows_ = 500;
const TString outDir_ = "plots";
// Selection
const double minJetE_ = 20;
const double minHt_ = 300.;
const double maxDeltaR_ = 0.2;


// === Global functions ===
void run();
void draw(const Data &data);
void drawResponsePdf(Response *dijetResp, Response *njetResp);
void fillResponse(const Data &data, Response *dijetResp, Response *njetResp);
Data readData(TChain *chain);
void selectData(Data &data);
void smearGenJets(Data &data, const Response *dijetResp, const Response *njetResp);


void dijetResolution() {
  run();
}


void run() {
  gErrorIgnoreLevel = 1001;
  gStyle->SetOptStat(0);

  TChain *chain = new TChain("DiJetTree");
  TString dir = "/scratch/hh/current/cms/user/stadie/Summer09QCDFlat-MC_31X_V9_7TeV-v1/";
  //TString dir = "~/lustre/mc/Summer09_QCDDiJet_AK5/Summer09_QCDDiJet_AK5_Pt0300-0380/";
  for(int i = 1; i <= 10; i++) {
    TString name = dir;
    name += "Summer09_QCD_Pt50_";
    //name += "Summer09_QCDDiJet_AK5_";
    name += i;
    name += ".root";
    chain->AddFile(name);
  }
  Data data = readData(chain);
  selectData(data);

  int nEBinEdges = 7;
  int nEtaBinEdges = 8;
  std::vector<double> eBinEdges(nEBinEdges);
  eBinEdges[0] = 10.;
  eBinEdges[1] = 40.;
  eBinEdges[2] = 80.;
  eBinEdges[3] = 120.;
  eBinEdges[4] = 200.;
  eBinEdges[5] = 500.;
  eBinEdges[6] = 2000.;

  std::vector<double> etaBinEdges(nEtaBinEdges);
  etaBinEdges[0] = -5.2;
  etaBinEdges[1] = -3.0;
  etaBinEdges[2] = -1.3;
  etaBinEdges[3] = -0.8;
  etaBinEdges[4] =  0.8;
  etaBinEdges[5] =  1.3;
  etaBinEdges[6] =  3.0;
  etaBinEdges[7] =  5.2;

  Response *dijetResp = new Response("DijetResponse",eBinEdges,etaBinEdges);
  Response *njetResp = new Response("NjetResponse",eBinEdges,etaBinEdges);
  fillResponse(data,dijetResp,njetResp);

  smearGenJets(data,dijetResp,njetResp);

  draw(data);  
  drawResponsePdf(dijetResp,njetResp);

  std::cout << "Binning of response:\n";
  for(int etaBin = 0; etaBin < nEtaBinEdges-1; etaBin++) {
    for(int eBin = 0; eBin < nEBinEdges-1; eBin++) {
      int idx = eBin + etaBin*(nEBinEdges-1);
      std::cout << " " << idx;
      std::cout << ":  E " << eBinEdges[eBin] << " - " << eBinEdges[eBin+1];
      std::cout << ",  eta " << etaBinEdges[etaBin] << " - " << etaBinEdges[etaBin+1] << std::endl;
    }
  }  
}


void draw(const Data &data) {
  std::cout << "Drawing histograms\n";
  
  // Pt spectra
  Hist hPtSpectrum(2);
  hPtSpectrum[0] = new TH1D("PtGenSpectrum",";p^{gen}_{T} (GeV)",50,0,2100);
  hPtSpectrum[1] = new TH1D("PtSpectrum",";p_{T} (GeV)",50,0,2100);

  TH1D *hMultiplicity = new TH1D("hMultiplicity",";Number of jets",20,0,20);

  std::vector<double> mhtBinEdges(15);
  for(int i = 0; i < 10; i++) {
    mhtBinEdges[i] = i*10.;
  }
  mhtBinEdges[10] = 120.;
  mhtBinEdges[11] = 150.;
  mhtBinEdges[12] = 200.;
  mhtBinEdges[13] = 300.;
  mhtBinEdges[14] = 600.;
  
  Hist hMHT(6);
  hMHT[0] = new TH1D("MHT",";MHT (GeV)",14,&(mhtBinEdges.front()));
  hMHT[1] = new TH1D("MHTGen",";Generator MHT (GeV)",14,&(mhtBinEdges.front()));
  hMHT[1]->SetLineColor(4);
  hMHT[2] = static_cast<TH1D*>(hMHT[0]->Clone("MHTPredDijet"));
  hMHT[2]->SetTitle("Prediction from dijets");
  hMHT[2]->SetLineColor(2);
  hMHT[3] = static_cast<TH1D*>(hMHT[0]->Clone("MHTPredNjet"));
  hMHT[3]->SetTitle("Prediction from njets");
  hMHT[3]->SetLineColor(4);
  hMHT[4] = static_cast<TH1D*>(hMHT[0]->Clone("MHTRatioDijet"));
  hMHT[4]->SetTitle("Prediction from dijets;MHT (GeV);Ratio prediction / simulation");
  hMHT[4]->SetLineColor(2);
  hMHT[5] = static_cast<TH1D*>(hMHT[0]->Clone("MHTRatioNjet"));
  hMHT[5]->SetTitle("Prediction from njets;MHT (GeV);Ratio prediction / simulation");
  hMHT[5]->SetLineColor(4);
  for(size_t i = 0; i < hMHT.size(); i++) {
    hMHT[i]->Sumw2();
  }    

  // Loop over data and fill hists
  for(DataIt it = data.begin(); it != data.end(); it++) {
    Event *evt = *it;
    
    if( evt->isSelected() ) {
      hMultiplicity->Fill(evt->nJets());

      for(JetIt jetIt = evt->firstJet(); jetIt != evt->lastJet(); jetIt++) {
	hPtSpectrum[0]->Fill((*jetIt)->ptGen());
	hPtSpectrum[1]->Fill((*jetIt)->pt());
      }

      hMHT[0]->Fill(evt->mht());
      hMHT[1]->Fill(evt->mhtGen());
      hMHT[2]->Fill(evt->mhtPredDijet());
      hMHT[3]->Fill(evt->mhtPredNjet());
      hMHT[4]->Fill(evt->mhtPredDijet());
      hMHT[5]->Fill(evt->mhtPredNjet());
     }
  } // End of loop over data
  hMHT[4]->Divide(hMHT[0]);
  hMHT[5]->Divide(hMHT[0]);

  TCanvas *can = new TCanvas("canMultiplicity","",500,500);
  can->cd();
  hMultiplicity->Draw();
  TString name = outDir_;
  name += "/JetMultiplicity.eps";
  can->SaveAs(name,"eps");

  for(size_t i = 0; i < hMHT.size(); i++) {
    can->cd();
    hMHT[i]->Draw();
    if( i < 4 ) can->SetLogy();
    name = outDir_;
    name += "/";
    name += hMHT[i]->GetName();
    name += ".eps";
    can->SaveAs(name,"eps");
    can->SetLogy(0);
  }

  can->cd();
  hMHT[4]->Draw();
  hMHT[5]->Draw("same");
  TLegend *leg = new TLegend(0.5,0.66,0.93,0.83);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  //  leg->SetTextSize(0.04);
  leg->AddEntry(hMHT[5],"Njet response","L");
  leg->AddEntry(hMHT[4],"Dijet response","L");
  leg->Draw("same");
  name = outDir_;
  name += "/MHTRatios.eps";
  can->SaveAs(name,"eps");

  Can canPtSpectrum(2);
  canPtSpectrum[0] = new TCanvas("canPtSpectrum0","GenJetPt",500,500);
  canPtSpectrum[1] = new TCanvas("canPtSpectrum1","JetPt",500,500);
  CanIt c = canPtSpectrum.begin();
  HistIt h = hPtSpectrum.begin();
  for(; c != canPtSpectrum.end(); c++, h++) {
    (*c)->cd();
    (*h)->SetNdivisions(505);
    (*h)->Draw();
    (*c)->SetLogy();
    name = outDir_;
    name += "/";
    name += (*h)->GetName();
    name += ".eps";
    (*c)->SaveAs(name,"eps");
  }
}


void drawResponsePdf(Response *dijetResp, Response *njetResp) {
  for(int c = 0; c < dijetResp->nBins(); c++) {
    TString name = "canDijetResp";
    name += c;
    TString title = "Dijet resp ";
    title += c;
    TCanvas *can = new TCanvas(name,title,500,500);
    can->cd();
    TH1D *h = dijetResp->pdf(c);
    h->GetYaxis()->SetRangeUser(0.1,10*h->GetMaximum());
    h->SetLineColor(2);
    h->Draw();
    can->SetLogy();
    name = outDir_;
    name += "/DijetResp_";
    name += c;
    name += ".eps";
    can->SaveAs(name,"eps");
    if( h->Integral() == 0 ) std::cout << "Dijet " << c << std::endl;
  }
  for(int c = 0; c < njetResp->nBins(); c++) {
    TString name = "canNjetResp";
    name += c;
    TString title = "Njet resp ";
    title += c;
    TCanvas *can = new TCanvas(name,title,500,500);
    can->cd();
    TH1D *h = njetResp->pdf(c);
    h->GetYaxis()->SetRangeUser(0.1,10*h->GetMaximum());
    h->SetLineColor(4);
    h->Draw();
    can->SetLogy();
    name = outDir_;
    name += "/NjetResp_";
    name += c;
    name += ".eps";
    can->SaveAs(name,"eps");
    if( h->Integral() == 0 ) std::cout << "Njet " << c << std::endl;
  }
  for(int c = 0; c < njetResp->nBins(); c++) {
    TString name = "canResp";
    name += c;
    TString title = "Resp ";
    title += c;
    TCanvas *can = new TCanvas(name,title,500,500);
    can->cd();
    TH1D *h1 = dijetResp->pdf(c);
    TH1D *h2 = njetResp->pdf(c);
    if( h1->Integral() ) h1->Scale(1./h1->Integral("width"));
    if( h2->Integral() ) h2->Scale(1./h2->Integral("width"));
    h1->GetYaxis()->SetRangeUser(2E-5,10*h1->GetMaximum());
    h1->Draw();
    h2->Draw("same");
    can->SetLogy();
    name = outDir_;
    name += "/Resp_";
    name += c;
    name += ".eps";
    can->SaveAs(name,"eps");
  }
}


void fillResponse(const Data &data, Response *dijetResp, Response *njetResp) {
  for(DataIt it = data.begin(); it != data.end(); it++) {
    Event *evt = *it;

    if( evt->isSelected() ) {
      int n = 0;
      for(JetIt jetIt = evt->firstJet(); jetIt != evt->lastJet(); jetIt++, n++) {
	Jet *jet = *jetIt;
	
 	if( njetResp->isInRange(jet->eGen(),jet->eta()) ) {
 	  njetResp->fill(jet->eGen(),jet->eta(),jet->e()/jet->eGen());
 	}
	if( evt->isDijet() && n < 2 && dijetResp->isInRange(jet->eGen(),jet->eta()) ) {
	  dijetResp->fill(jet->eGen(),jet->eta(),jet->e()/jet->eGen());
	}
      }
    }
  }
}


Data readData(TChain *chain) {
  std::cout << "Reading data\n";;

  Data data;
    
  // Init read quantities
  int nObjJet = 0;
  int nObjGenJet = 0;
  int nObjTow = 0;             
  float genJetColE[maxNJet_];
  float genJetColPt[maxNJet_];
  float genJetColEta[maxNJet_];
  float genJetColPhi[maxNJet_];
  float jetE[maxNJet_];
  float jetPt[maxNJet_];
  float jetEta[maxNJet_];
  float jetPhi[maxNJet_];
  float jetCorrL2L3[maxNJet_];
  int genJetColJetIdx[maxNJet_];
  int towJetIdx[maxNTows_];
  float towEm[maxNTows_];

  // Set branch addresses
  chain->SetBranchAddress("NobjJet",&nObjJet);
  chain->SetBranchAddress("NobjGenJet",&nObjGenJet);
  chain->SetBranchAddress("NobjTow",&nObjTow);
  chain->SetBranchAddress("GenJetColE",genJetColE);
  chain->SetBranchAddress("GenJetColPt",genJetColPt);
  chain->SetBranchAddress("GenJetColEta",genJetColEta);
  chain->SetBranchAddress("GenJetColPhi",genJetColPhi);
  chain->SetBranchAddress("JetE",jetE);
  chain->SetBranchAddress("JetPt",jetPt);
  chain->SetBranchAddress("JetEta",jetEta);
  chain->SetBranchAddress("JetPhi",jetPhi);
  chain->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  chain->SetBranchAddress("GenJetColJetIdx",genJetColJetIdx);
  chain->SetBranchAddress("Tow_jetidx",towJetIdx);
  chain->SetBranchAddress("TowEm",towEm);
  
  // Loop over tree entries
  for(int n = 0; n < chain->GetEntries(); n++) {
    chain->GetEntry(n);
    
    if( nObjJet > maxNJet_ || nObjGenJet > maxNJet_ ) {
      std::cerr << "WARNING: nObjJet = " << nObjJet << " > maxNJet_. Skipping event.\n";
      continue;
    }
    if( nObjTow > maxNTows_ ) {
      std::cerr << "WARNING: nObjTow = " << nObjTow << " > maxNTows_. Skipping event.\n";
      continue;
    }

    // Create Event and fill jets
    Event *evt = new Event();
    // Loop over gen jets
    for(int genJetIdx = 0; genJetIdx < nObjGenJet; genJetIdx++) {
      int jetIdx = genJetColJetIdx[genJetIdx];

      double e = jetCorrL2L3[jetIdx]*jetE[jetIdx];
      double pt = jetCorrL2L3[jetIdx]*jetPt[jetIdx];
      double ptGen = genJetColPt[genJetIdx];
      double eGen = genJetColE[genJetIdx];
      if( eGen < minJetE_ ) continue;

      // Calculate emf from towers in this jet
      double emf = 0.;
      for(int t = 0; t < nObjTow; t++)
        {
          if( towJetIdx[t] == jetIdx ) {
	    emf += towEm[t];
	  }
	}
      emf /= jetPt[jetIdx];
      
      // Calculate dR between genjet and jet
      double dEta = genJetColEta[genJetIdx] - jetEta[jetIdx];
      double dPhi = genJetColPhi[genJetIdx] - jetPhi[jetIdx];
      double dR = sqrt(dEta*dEta + dPhi*dPhi);

      // Add jet to event
      evt->addJet(new Jet(eGen,ptGen,e,pt,jetEta[jetIdx],jetPhi[jetIdx],emf,dR));      
    } // End of loop over gen jets
    data.push_back(evt);
  } // End of loop over tree entries
  std::cout << " Read " << data.size() << " events\n";

  return data;
}


void selectData(Data &data) {
  std::cout << "Selecting events\n";
  int nSelectedEvts = 0;
  int nDijetEvts = 0;
  std::vector<Event*>::iterator it = data.begin();
  for(; it != data.end(); it++) {
    Event *evt = *it;

    //    if( evt->htGen() < minHt_ ) continue;
    if( evt->maxDeltaR() > maxDeltaR_ ) continue;
    nSelectedEvts++;
    evt->isSelected(true);

    if( evt->rel3rdJetPt() > 0.1 || std::abs(evt->deltaPhi()) < 2.7 ) continue;
    nDijetEvts++;
    evt->isDijet(true);
  }
  std::cout << " Selected " << nSelectedEvts << " events\n";
  std::cout << " Selected " << nDijetEvts << " dijet events\n";
}


void smearGenJets(Data &data, const Response *dijetResp, const Response *njetResp) {
  std::vector<Event*>::iterator it = data.begin();
  for(; it != data.end(); it++) {
    Event *evt = *it;

    bool isSelected = true;
    double mhtDijet[2] = { 0., 0. };
    double mhtNjet[2] = { 0., 0. };
    for(JetIt jetIt = evt->firstJet(); jetIt != evt->lastJet(); jetIt++) {
      Jet *jet = *jetIt;

      double ptSmearDijet = dijetResp->random(jet->eGen(),jet->eta())*jet->ptGen();
      double ptSmearNjet = njetResp->random(jet->eGen(),jet->eta())*jet->ptGen();
      if( ptSmearDijet < 0 || ptSmearNjet < 0 ) isSelected = false;

      mhtDijet[0] += ptSmearDijet*cos(jet->phi());
      mhtDijet[1] += ptSmearDijet*sin(jet->phi());
      mhtNjet[0] += ptSmearNjet*cos(jet->phi());
      mhtNjet[1] += ptSmearNjet*sin(jet->phi());
    }
    evt->mhtPredDijet(sqrt( mhtDijet[0]*mhtDijet[0] + mhtDijet[1]*mhtDijet[1] ));
    evt->mhtPredNjet(sqrt( mhtNjet[0]*mhtNjet[0] + mhtNjet[1]*mhtNjet[1] ));  
    if( !isSelected ) evt->isSelected(false);
  }
}
