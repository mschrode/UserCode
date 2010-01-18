// $Id: firstDataAnalysis.C,v 1.8 2010/01/17 15:21:41 mschrode Exp $

// ===== Script for dijet analysis from Kalibri ntuples =====
//
// This script takes Kalibri ntuples as an input, performes
// a dijet selection, and plots some quantities of the selected
// events. The n-1 distributions of the selection are also
// plotted and a cut flow is printed (optionally into a LaTeX
// table). The distributions of different samples can be 
// superimposed (normalized to the first sample) in order to
// e.g. compare MC with data.
//
// Author: Matthias Schroeder
// Date: So 17 Jan 2010 16:24:25 CET

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TVector2.h"



// ===== Type declarations =====
// =============================

class Jet;
class Event;
class Cut;
class CutFlow;
class Sample;

typedef std::vector<Event*> Data;
typedef std::vector<Event*>::iterator DataIt;

typedef std::vector<TH1D*> Distributions;
typedef std::vector<TH1D*>::iterator DistIt;




// ===== Class implementations =====
// =================================

// ----- Jet -----
class Jet {
 public:
  static bool ptGreaterThan(const Jet *j1, const Jet *j2);
  static bool corrPtGreaterThan(const Jet *j1, const Jet *j2);

  Jet(double pt,
      double eta,
      double phi,
      double emf,
      int n90Hits,
      double fHPD,
      double fRBX,
      double corrL2L3);

  double pt() const { return pt_; }
  double eta() const { return eta_; }
  double phi() const { return phi_; }
  double emf() const { return emf_; }
  int n90Hits() const { return n90Hits_; }
  double fHPD() const { return fHPD_; }
  double fRBX() const { return fRBX_; }
  double corrL2L3Pt() const { return corrL2L3Pt_; }

 private:
  const double pt_;
  const double eta_;
  const double phi_;
  const double emf_;
  const int n90Hits_;
  const double fHPD_;
  const double fRBX_;
  const double corrL2L3Pt_;
};

Jet::Jet(double pt,
	 double eta,
	 double phi,
	 double emf,
	 int n90Hits,
	 double fHPD,
	 double fRBX,
	 double corrL2L3)
  : pt_(pt),
    eta_(eta),
    phi_(phi),
    emf_(emf),
    n90Hits_(n90Hits),
    fHPD_(fHPD),
    fRBX_(fRBX),
    corrL2L3Pt_(corrL2L3*pt) {};

//! For sorting jets in calo pt
bool Jet::ptGreaterThan(const Jet *j1, const Jet *j2) {
  // check for 0
  if (j1 == 0) {
    return j2 != 0;
  } else if (j2 == 0) {
    return false;
  } else {
    return j1->pt() > j2->pt();
  }
}

//! For sorting jets in corrected pt
bool Jet::corrPtGreaterThan(const Jet *j1, const Jet *j2) {
  // check for 0
  if (j1 == 0) {
    return j2 != 0;
  } else if (j2 == 0) {
    return false;
  } else {
    return j1->corrL2L3Pt() > j2->corrL2L3Pt();
  }
}


// ----- Event -----
class Event {
 public:
  Event(unsigned int runNumber,
	unsigned int lumiBlockNumber,
	unsigned int vtxNTracks,
	double vtxPosZ,
	double sumEt,
	double met,
	const std::vector<Jet*>& jets);
  ~Event();

  unsigned int runNumber() const { return runNumber_; }
  unsigned int lumiBlockNumber() const { return lumiBlockNumber_; }
  unsigned int vtxNTracks() const { return vtxNTracks_; }
  double vtxPosZ() const { return vtxPosZ_; }
  double sumEt() const { return sumEt_; }
  double met() const { return met_; }
  int nJets() const { return static_cast<int>(jets_.size()); }
  const Jet *jet(int idx) const { return idx >= 0 && idx < nJets() ? jets_[idx] : 0; }
  const Jet *corrJet(int idx) const { return idx >= 0 && idx < nJets() ? corrJets_[idx] : 0; }

 private:
  const unsigned int runNumber_;
  const unsigned int lumiBlockNumber_;
  const unsigned int vtxNTracks_;
  const double vtxPosZ_;
  const double sumEt_;
  const double met_;
  std::vector<Jet*> jets_;
  std::vector<Jet*> corrJets_;
};

Event::Event(unsigned int runNumber,
	     unsigned int lumiBlockNumber,
	     unsigned int vtxNTracks,
	     double vtxPosZ,
	     double sumEt,
	     double met,
	     const std::vector<Jet*>& jets)
  : runNumber_(runNumber),
    lumiBlockNumber_(lumiBlockNumber),
    vtxNTracks_(vtxNTracks),
    vtxPosZ_(vtxPosZ),
    sumEt_(sumEt),
    met_(met),
    jets_(jets),
    corrJets_(jets) {
  // Sort jets by pt
  std::sort(jets_.begin(),jets_.end(),Jet::ptGreaterThan);
  // Sort corrected jets by L2L3 corrected pt
  std::sort(corrJets_.begin(),corrJets_.end(),Jet::corrPtGreaterThan);  
}

Event::~Event() {
  std::vector<Jet*>::iterator jetIt = jets_.begin();
  for(; jetIt != jets_.end(); jetIt++) {
    delete *jetIt;
  }
  jets_.clear();
  corrJets_.clear();
}


// ----- Cut -----
class Cut {
public:
  Cut(const TString &sampleName, int idx, const TString &id, const TString &label, bool isEvtCut, double min, double max, bool abs, int nBins, double xMin, double xMax, const TString &xAxisTitle,bool logy = false);
  Cut(const TString &sampleName, int idx, const TString &id, const TString &label, bool isEvtCut, double val, bool low, bool abs, int nBins, double xMin, double xMax, const TString &xAxisTitle, bool logy = false);
  ~Cut();

  bool abs() const { return abs_; }
  bool isEvtCut() const { return isEvtCut_; }
  bool hasMin() const { return hasMin_; }
  bool hasMax() const { return hasMax_; }
  double min() const { return min_; }
  double max() const { return max_; }
  int nPassed() const { return nPassed_; }
  int nRejected() const { return nRejected_; }
  TString distName() const;
  bool distLogY() const { return logY_; }
  TString printLineScreen() const { return printLine(false); }
  TString printLineLaTeX() const { return printLine(true); }

  TH1D *dist() { return dist_; }
  bool passes(double val1, double val2 = 0.);

private:
  const int idx_;
  const TString id_;
  const TString label_;
  const bool isEvtCut_; // Cut on event or on first two jets
  const double min_;
  const double max_;
  const bool abs_;
  const bool hasMin_;
  const bool hasMax_;
  const bool logY_;

  TH1D *dist_;
  int nPassed_;
  int nRejected_;

  void init(const TString &sampleName, int nBins, double xMin, double xMax, const TString &xAxisTitle);
  TString printLine(bool latex) const;
};

Cut::Cut(const TString &sampleName, int idx, const TString &id, const TString &label, bool isEvtCut, double min, double max, bool abs, int nBins, double xMin, double xMax, const TString &xAxisTitle, bool logy)
  : idx_(idx), id_(id), label_(label), isEvtCut_(isEvtCut), min_(min), max_(max), abs_(abs), hasMin_ (true), hasMax_(true), logY_(logy) {
  init(sampleName,nBins,xMin,xMax,xAxisTitle);
}

Cut::Cut(const TString &sampleName, int idx, const TString &id, const TString &label, bool isEvtCut, double val, bool low, bool abs, int nBins, double xMin, double xMax, const TString &xAxisTitle, bool logy)
  : idx_(idx), id_(id), label_(label), isEvtCut_(isEvtCut), min_(val), max_(val), abs_(abs), hasMin_ (low), hasMax_(!low), logY_(logy) {
  init(sampleName,nBins,xMin,xMax,xAxisTitle);
}

Cut::~Cut() {
  delete dist_;
}

bool Cut::passes(double val1, double val2) {
  bool passes = true;
  int nVal = 2;
  if( isEvtCut() ) nVal = 1;
  for(int i = 0; i < nVal; i++) {
    bool iPass = false;
    double val = val1;
    if( i == 1 ) val = val2;
    dist_->Fill(val);
    if( abs() ) val = std::abs(val);

    if( hasMin() && hasMax() ) {
      if( val >= min() && val <= max() ) iPass = true;
    } else if( hasMin() ) {
      if( val >= min() ) iPass = true;
    } else if( hasMax() ) {
      if( val <= max() ) iPass = true;
    }
    passes = passes && iPass;
  }

  if( passes ) nPassed_++;
  else nRejected_++;

  return passes;
}

TString Cut::distName() const {
  TString name = "Cut_";
  name += idx_;
  name += "_";
  name += id_;
  return name;
}

TString Cut::printLine(bool latex) const {
  TString line;
  if( idx_ < 10 ) line += " ";
  line += idx_;
  if( latex ) line += " & \\texttt{";
  else line += ": ";
  if( abs() ) line += "|";
  line += label_;
  if( abs() ) line += "|";
  if( latex ) line += "}";
  if( hasMin() ) {
    char min[20];
    sprintf(min,"%.2f",min_);
    if( latex ) line += " $ \\geq ";
    else line += " >= ";
    line += min;
  }
  if( hasMax() ) {
    char max[20];
    sprintf(max,"%.2f",max_);
    if( latex ) line += " $ \\leq ";
    else line += " <= ";
    line += max;
  }
  // In case of int remove decimal digits
  size_t pos = line.First('.');
  TString end = line(pos+1,line.Length()-1);
  while( end.Length() > 0 ) {
    if( end.Chop() != "0" ) break;
  }
  if( end.Length() == 0 ) {
    line = line(0,pos);
  }
  if( latex ) line += " $ ";
  else while( line.Length() < 40 ) line += " ";

  return line;
}

void Cut::init(const TString &sampleName, int nBins, double xMin, double xMax, const TString &xAxisTitle) {
  TString name = sampleName;
  name += ":Cut_";
  name += idx_;
  name += "_";
  name += id_;
  TString title = "Cut ";
  title += idx_;
  title += ": ";
  title += label_;
  title += " (n-1 plot)";
  dist_ = new TH1D(name,title,nBins,xMin,xMax);
  dist_->SetXTitle(xAxisTitle);
  if( isEvtCut() ) dist_->SetYTitle("Number of events");
  else dist_->SetYTitle("Number of jets");
  dist_->SetMarkerStyle(20);
  
  nPassed_ = 0;
  nRejected_ = 0;
}

// ----- CutFlow -----
class CutFlow {
public:
  CutFlow(const TString &sampleName);
  ~CutFlow();

  bool abs(int idx) const { return isValidCutIdx(idx) ? cuts_[idx]->abs() : false; }
  double min(int idx) const { return isValidCutIdx(idx) ? cuts_[idx]->min() : 0; }
  double max(int idx) const { return isValidCutIdx(idx) ? cuts_[idx]->max() : 0; }
  double hasMin(int idx) const { return isValidCutIdx(idx) ? cuts_[idx]->hasMin() : 0; }
  double hasMax(int idx) const { return isValidCutIdx(idx) ? cuts_[idx]->hasMax() : 0; }
  int nTotalEvts() const { 
    return isValidCutIdx(0) ? cuts_[0]->nRejected() + cuts_[0]->nPassed() : 0;
  }
  int nPassedEvts() const { 
    return isValidCutIdx(nCuts()-1) ? cuts_[nCuts()-1]->nPassed() : 0;
  }
  int nPassedEvts(int idx) const { 
    return isValidCutIdx(idx) ? cuts_[idx]->nPassed() : 0;
  }
  int nCuts() const { return static_cast<int>(cuts_.size()); }
  TString distName(int idx) const { return isValidCutIdx(idx) ? cuts_[idx]->distName() : ""; }
  bool distLogY(int idx) const { return isValidCutIdx(idx) ? cuts_[idx]->distLogY() : false; }
  double distIntegral(int idx) const { 
    return isValidCutIdx(idx) ? cuts_[idx]->dist()->Integral() : 0;
  }
  TString printLineScreen(int idx) const {
    return isValidCutIdx(idx) ? cuts_[idx]->printLineScreen() : "";
  }
  TString printLineLaTeX(int idx) const {
    return isValidCutIdx(idx) ? cuts_[idx]->printLineLaTeX() : "";
  }
  
  TH1D *dist(int idx) { return isValidCutIdx(idx) ? cuts_[idx]->dist() : 0; }
  void normaliseDists(const CutFlow *cutFlow);
  bool passes(const Event *evt);

private:
  std::vector<Cut*> cuts_;

  bool isValidCutIdx(int idx) const { return (idx >= 0 && idx < nCuts() ); }
};

CutFlow::CutFlow(const TString &sampleName) {
  cuts_.push_back(new Cut(sampleName,0,"NJets","N jets",true,2.,true,false,
			  30,0.,30.,"Number of jets",true));

  cuts_.push_back(new Cut(sampleName,1,"NVtxTracks","N vtx tracks",true,2.,true,false,
			  50,0.,50.,"Number of primary vertex tracks"));

  cuts_.push_back(new Cut(sampleName,2,"VtxZPos","Vtx z",true,20.,false,true,
			  30,-30.,30.,"Primary vertex z (cm)"));

  cuts_.push_back(new Cut(sampleName,3,"RelMET","MET / SumEt",true,0.5,false,false,
			  50,0.,1.,"MET / sumEt",true));

  cuts_.push_back(new Cut(sampleName,4,"PtJet2","Pt jet (2)",true,4.,true,false,
			  30,0.,30.,"p^{jet2}_{T} (GeV)",true));

  cuts_.push_back(new Cut(sampleName,5,"Eta","Eta jet (1,2)",false,3.,false,true,
			  20,-5.2,5.2,"#eta"));

  cuts_.push_back(new Cut(sampleName,6,"DeltaPhi","Delta phi jet (1,2)",true,2.1,true,true,
			  20,-3.2,3.2,"#Delta#phi(jet1,jet2)"));

  cuts_.push_back(new Cut(sampleName,7,"fEMF","fEMF jet (1,2)",false,0.01,true,false,
			  20,0.,1.,"f_{EM}"));

  cuts_.push_back(new Cut(sampleName,8,"N90Hits","n90Hits jet (1,2)",false,2,true,false,
			  20,0.,50.,"n90Hits"));

  cuts_.push_back(new Cut(sampleName,9,"fHPD","fHPDX jet (1,2)",false,0.98,false,false,
			  20,0.,1.,"f_{HPD}"));

  cuts_.push_back(new Cut(sampleName,10,"fRBX","fRBX jet (1,2)",false,0.98,false,false,
			  20,0.,1.,"f_{RBX}"));

}

CutFlow::~CutFlow() {
  for(std::vector<Cut*>::iterator it = cuts_.begin();
      it != cuts_.end(); it++) {
    delete *it;
  }
  cuts_.clear();
}

void CutFlow::normaliseDists(const CutFlow *cutFlow) {
  for(int i = 0; i < nCuts(); i++) {
    double norm = cutFlow->distIntegral(i);
    if( distIntegral(i) ) {
      norm /= distIntegral(i);
      cuts_[i]->dist()->Scale(norm);
    }
  }
}

bool CutFlow::passes(const Event *evt) {
  bool passes = true;
  if( !cuts_[0]->passes( evt->nJets()) ) passes = false;
  else if( !cuts_[1]->passes(evt->vtxNTracks()) ) passes = false;
  else if( !cuts_[2]->passes(evt->vtxPosZ()) ) passes = false;
  else if( !cuts_[3]->passes(evt->met()/evt->sumEt()) ) passes = false;
  else if( !cuts_[4]->passes(evt->jet(1)->pt()) ) passes = false;
  else if( !cuts_[5]->passes(evt->jet(0)->eta(),evt->jet(1)->eta()) ) passes = false;
  else if( !cuts_[6]->passes(TVector2::Phi_mpi_pi(evt->jet(0)->phi()-evt->jet(1)->phi())) ) passes = false;
  else if( !cuts_[7]->passes(evt->jet(0)->emf(),evt->jet(1)->emf()) ) passes = false;
  else if( !cuts_[8]->passes(evt->jet(0)->n90Hits(),evt->jet(1)->n90Hits()) ) passes = false;
  else if( !cuts_[9]->passes(evt->jet(0)->fHPD(),evt->jet(1)->fHPD()) ) passes = false;
  else if( !cuts_[10]->passes(evt->jet(0)->fRBX(),evt->jet(1)->fRBX()) ) passes = false;

  return passes;
}


// ----- Sample -----
class Sample {
public:
  Sample(const TString &sampleName, const TString &treeName);
  ~Sample();

  TString name() const { return name_; }
  TString drawOptions() const { return drawOpt_; }

  int nDists() const { return static_cast<int>(dists_.size()); }
  TString distName(int idx) const;
  bool distLogY(int idx) const { return isValidDistIdx(idx) ? distLogY_[idx] : false; }
  void normaliseDists(Sample *sample);
  double distIntegral(int idx) const { return isValidDistIdx(idx) ? dists_[idx]->Integral() : 0; }
  TH1D *dist(int idx) const { return isValidDistIdx(idx) ? dists_[idx]: 0; }

  CutFlow *cutFlow() { return cutFlow_; }
  int nTotalEvts() const { return cutFlow_->nTotalEvts(); }
  int nPassedEvts() const { return cutFlow_->nPassedEvts(); }
  int nPassedEvts(int idx) const { return cutFlow_->nPassedEvts(idx); }
  int nCuts() const { return cutFlow_->nCuts(); }
  TH1D *cutDist(int idx) { return cutFlow_->dist(idx); }
  TString cutDistName(int idx) const { return cutFlow_->distName(idx); }
  bool cutDistLogY(int idx) const { return cutFlow_->distLogY(idx); }
  TString printCutLineScreen(int idx) const { return cutFlow_->printLineScreen(idx); }
  TString printCutLineLaTeX(int idx) const { return cutFlow_->printLineLaTeX(idx); }

  void addFile(const TString &name);
  void fillDistributions();
  void readData(int nMaxEvts = -1);
  void setDrawOptions(const TString &opt) { drawOpt_ = opt; }

private:
  const TString name_;
  const int maxNJet_;

  TChain *chain_;
  Data data_;
  Distributions dists_;
  std::vector<bool> distLogY_;
  std::vector<unsigned int> runs_;
  TString drawOpt_;
  CutFlow *cutFlow_;
  TRandom3 *rand_;

  bool isValidDistIdx(int idx) const { return (idx >= 0 && idx < nDists() ); }
};

Sample::Sample(const TString &sampleName,const TString &treeName)
  : name_(sampleName), maxNJet_(50) {
  chain_ = new TChain(treeName);
  dists_ = Distributions(6);
  distLogY_ = std::vector<bool>(dists_.size(),false);

  TString name = name_;
  name += ":Pt";
  dists_[0] = new TH1D(name,"Selected dijet events",15,0,60);
  dists_[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
  dists_[0]->GetYaxis()->SetTitle("Number of jets");
  distLogY_[0] = true;

  name = name_;
  name += ":PtCorr";
  dists_[1] = new TH1D(name,"Selected dijet events",15,0,60);
  dists_[1]->GetXaxis()->SetTitle("Corrected p_{T} (GeV)");
  dists_[1]->GetYaxis()->SetTitle("Number of jets");
  distLogY_[1] = true;

  name = name_;
  name += ":PtAsym";
  dists_[2] = new TH1D(name,"Selected dijet events",15,-1.5,1.5);
  dists_[2]->GetXaxis()->SetTitle("p_{T} asymmetry");
  dists_[2]->GetYaxis()->SetTitle("Number of events");

  name = name_;
  name += ":Phi";
  dists_[3] = new TH1D(name,"Selected dijet events",15,-3.,3.);
  dists_[3]->GetXaxis()->SetTitle("#phi");
  dists_[3]->GetYaxis()->SetTitle("Number of jets");

  name = name_;
  name += ":DeltaPhi";
  dists_[4] = new TH1D(name,"Selected dijet events",15,2.,4.);
  dists_[4]->GetXaxis()->SetTitle("#Delta#phi");
  dists_[4]->GetYaxis()->SetTitle("Number of events");

  name = name_;
  name += ":Eta";
  dists_[5] = new TH1D(name,"Selected dijet events",15,-5,5);
  dists_[5]->GetXaxis()->SetTitle("#eta");
  dists_[5]->GetYaxis()->SetTitle("Number of jets");

  for(int d = 0; d < nDists(); d++) {
    dists_[d]->SetMarkerStyle(20);
  }

  cutFlow_ = new CutFlow(name_);
  rand_ = new TRandom3(0);
}

Sample::~Sample() {
  delete chain_;
  for(DataIt it = data_.begin(); it != data_.end(); it++) {
    delete *it;
  }
  data_.clear();
  for(DistIt it = dists_.begin(); it != dists_.end(); it++) {
    delete *it;
  }
  dists_.clear();
  runs_.clear();
  delete cutFlow_;
  delete rand_;
}

TString Sample::distName(int idx) const {
  TString name;
  if( isValidDistIdx(idx) ) {
    name = dists_[idx]->GetName();
    size_t pos = name.First(':');
    name = name(pos+1,name.Length()-pos);
  }
  return name;
}

void Sample::addFile(const TString &name) {
  std::cout << "Sample '" << name_ << "': Adding file '" << name << "... " << std::flush;
  chain_->AddFile(name);
  std::cout << "ok\n";
}

void Sample::fillDistributions() {
  std::cout << "Sample '" << name_ << "': Filling distributions... " << std::flush;

  // Reset distributions
  for(DistIt it = dists_.begin(); it != dists_.end(); it++) {
    (*it)->Reset();
  }

  // Loop over data
  for(DataIt it = data_.begin(); it != data_.end(); it++) {
    Event *evt = *it;

    double diff = evt->jet(0)->pt() - evt->jet(1)->pt();
    if( rand_->Uniform() > 0.5 ) diff = evt->jet(1)->pt() - evt->jet(0)->pt();
    dists_[2]->Fill(diff/(evt->jet(0)->pt() + evt->jet(1)->pt()));

    double deltaPhi = TVector2::Phi_mpi_pi(evt->jet(0)->phi() - evt->jet(1)->phi());
    if( deltaPhi < 0 ) deltaPhi += 2*M_PI;
    dists_[4]->Fill(deltaPhi);

    // Loop over two jets leading in pt
    for(int j = 0; j < 2; j++) {
      dists_[0]->Fill(evt->jet(j)->pt());
      dists_[1]->Fill(evt->jet(j)->corrL2L3Pt());
      dists_[3]->Fill(evt->jet(j)->phi());
      dists_[5]->Fill(evt->jet(j)->eta());
    } // End of loop over dijets
  } // End of loop over data

  std::cout << "ok" << std::endl;
}

void Sample::normaliseDists(Sample *sample) {
  for(int i = 0; i < nDists(); i++) {
    double norm = sample->distIntegral(i);
    if( distIntegral(i) ) {
      norm /= distIntegral(i);
      dists_[i]->Scale(norm);
    }
  }
  cutFlow_->normaliseDists(sample->cutFlow());
}

void Sample::readData(int nMaxEvts) {
  std::cout << "Sample '" << name_ << "': Reading data... " << std::flush;

    // Reset data
    for(DataIt it = data_.begin(); it != data_.end(); it++) {
      delete *it;
    }
    data_.clear();

    // Init read quantities
    unsigned int runNumber = 0;
    unsigned int lumiBlockNumber = 0;
    int vtxNTracks = 0;
    float vtxPosZ = 0.;
    float sumEt = 0.;
    float met = 0.;
    int nObjJet = 0;
    float jetPt[maxNJet_];
    float jetEta[maxNJet_];
    float jetPhi[maxNJet_];
    float jetEMF[maxNJet_];
    int jetN90Hits[maxNJet_];
    float jetFHPD[maxNJet_];
    float jetFRBX[maxNJet_];
    float jetCorrL2L3[maxNJet_];

    // Set branch addresses
    chain_->SetBranchAddress("RunNumber",&runNumber);
    chain_->SetBranchAddress("LumiBlockNumber",&lumiBlockNumber);
    chain_->SetBranchAddress("VtxNTracks",&vtxNTracks);
    chain_->SetBranchAddress("VtxPosZ",&vtxPosZ);
    chain_->SetBranchAddress("Met",&met);
    chain_->SetBranchAddress("MetSum",&sumEt);
    chain_->SetBranchAddress("NobjJet",&nObjJet);
    chain_->SetBranchAddress("JetPt",jetPt);
    chain_->SetBranchAddress("JetEta",jetEta);
    chain_->SetBranchAddress("JetPhi",jetPhi);
    chain_->SetBranchAddress("JetEMF",jetEMF);
    chain_->SetBranchAddress("JetN90Hits",jetN90Hits);
    chain_->SetBranchAddress("JetFHPD",jetFHPD);
    chain_->SetBranchAddress("JetFRBX",jetFRBX);
    chain_->SetBranchAddress("JetCorrL2L3",jetCorrL2L3);
  
    // Loop over tree entries
    int nEntries = chain_->GetEntries();
    if( nMaxEvts > 0 && nEntries > nMaxEvts ) nEntries = nMaxEvts;
    for(int n = 0; n < nEntries; n++) {
      chain_->GetEntry(n);

      if( nObjJet > maxNJet_ ) {
	std::cerr << "WARNING: nObjJet = " << nObjJet << " > maxNJet_. Skipping event.\n";
	continue;
      }

      //      if( nObjJet < 2 ) continue;

      // Create event
      std::vector<Jet*> jets(nObjJet);
      for(int i = 0; i < nObjJet; i++) {
	jets[i] = new Jet(jetPt[i],jetEta[i],jetPhi[i],jetEMF[i],jetN90Hits[i],
			  jetFHPD[i],jetFRBX[i],jetCorrL2L3[i]);
      }
      Event *evt = new Event(runNumber,lumiBlockNumber,vtxNTracks,vtxPosZ,sumEt,met,jets);
      dummy++;
      if( cutFlow_->passes(evt) ) {
	data_.push_back(evt);
      }

      // Run number
      bool isNewRun = true;
      for(std::vector<unsigned int>::const_iterator r = runs_.begin();
	  r != runs_.end(); r++) {
	if( runNumber == *r ) {
	  isNewRun = false;
	  break;
	}
      }
      if( isNewRun ) {
	runs_.push_back(runNumber);
      }
    } // End of loop over tree entries
    std::sort(runs_.begin(),runs_.end());
}




// ===== Main script =====
// =======================

void draw(std::vector<Sample*> samples, const TString &label, const TString &prefix) {
  if( samples.size() > 0 ) {
    std::cout << "Drawing distributions... " << std::flush;
    gStyle->SetOptStat(0);

    // Normalize distributions to sample 0
    for(size_t i = 1; i < samples.size(); i++) {
      samples[i]->normaliseDists(samples[0]);
    }

    // Draw and save cut flow distributions
    for(int d = 0; d < samples[0]->nCuts(); d++) {
      TString canName = "can";
      canName += samples[0]->cutDistName(d);
      TCanvas *can = new TCanvas(canName,samples[0]->cutDistName(d),500,500);
      can->cd();

      TLegend *leg = new TLegend(0.5,0.85-samples.size()*0.08,0.93,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(10);
      leg->SetTextFont(42);
      leg->SetHeader(label);

      std::vector<TLine*> lines;

      for(size_t s = 0; s < samples.size(); s++) {
	TString opt = samples[s]->drawOptions();
	if( s > 0 ) opt += "same";
	TH1D *h = samples[s]->cutDist(d);
	if( h ) h->Draw(opt);
	if( s == 0 ) {
	  char text[50];
      	  sprintf(text,"%s: %i entries",samples[s]->name().Data(),static_cast<int>(h->GetEntries()));
	  leg->AddEntry(h,text,"PL");

	  double min = 0.;
	  double max = 2.*h->GetMaximum();
	  if( samples[s]->cutDistLogY(d) ) {
	    min = 0.1;
	    max *= 2.5;
	  }
	  h->GetYaxis()->SetRangeUser(min,max);

	  if( samples[0]->cutFlow()->hasMin(d) ) {
	    lines.push_back(new TLine(samples[0]->cutFlow()->min(d),min,
				      samples[0]->cutFlow()->min(d),max));
	    if( samples[0]->cutFlow()->abs(d) ) {
	      lines.push_back(new TLine(-samples[0]->cutFlow()->min(d),min,
					-samples[0]->cutFlow()->min(d),max));
	    }
	  }
	  if( samples[0]->cutFlow()->hasMax(d) ) {
	    lines.push_back(new TLine(samples[0]->cutFlow()->max(d),min,
				      samples[0]->cutFlow()->max(d),max));
	    if( samples[0]->cutFlow()->abs(d) ) {
	      lines.push_back(new TLine(-samples[0]->cutFlow()->max(d),min,
					-samples[0]->cutFlow()->max(d),max));
	    }
	  }
	} else {
	  char text[50];
      	  sprintf(text,"%s: normalized",samples[s]->name().Data());
	  leg->AddEntry(h,"MC: normalized","L");
	}
      }
      leg->Draw("same");
      if( samples[0]->cutDistLogY(d) ) can->SetLogy();

      for(std::vector<TLine*>::iterator it = lines.begin();
	  it != lines.end(); it++) {
	(*it)->SetLineWidth(2);
	(*it)->SetLineColor(2);
	(*it)->SetLineStyle(2);
	(*it)->Draw("same");
      }

      TString fileName = prefix;
      fileName += "_";
      fileName += samples[0]->cutDistName(d);
      fileName += ".eps";
      can->SaveAs(fileName,"eps");
    }

    // Draw and save distributions
    for(int d = 0; d < samples[0]->nDists(); d++) {
      TString canName = "can";
      canName += samples[0]->distName(d);
      TCanvas *can = new TCanvas(canName,samples[0]->distName(d),500,500);
      can->cd();

      TLegend *leg = new TLegend(0.5,0.85-samples.size()*0.08,0.93,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetHeader(label);

      for(size_t s = 0; s < samples.size(); s++) {
	TString opt = samples[s]->drawOptions();
	if( s > 0 ) opt += "same";
	TH1D *h = samples[s]->dist(d);
	if( h ) h->Draw(opt);
	if( s == 0 ) {
	  char text[50];
      	  sprintf(text,"%s: %i entries",samples[s]->name().Data(),static_cast<int>(h->GetEntries()));
	  leg->AddEntry(h,text,"PL");

	  double max = h->GetMaximum();
	  h->GetYaxis()->SetRangeUser(0.,2.*max);
	  if( samples[s]->distLogY(d) ) {
	    h->GetYaxis()->SetRangeUser(0.1,5*h->GetMaximum());
	  }
	} else {
	  char text[50];
      	  sprintf(text,"%s: normalized",samples[s]->name().Data());
	  leg->AddEntry(h,"MC: normalized","L");
	}
      }
      leg->Draw("same");
      if( samples[0]->distLogY(d) ) can->SetLogy();

      TString fileName = prefix;
      fileName += "_";
      fileName += samples[0]->distName(d);
      fileName += ".eps";
      can->SaveAs(fileName,"eps");
    }
    std::cout << "ok" << std::endl;
  }
}

void printCutFlow(const std::vector<Sample*> samples) {
  if( samples.size() > 0 ) {
    std::cout << "\n\n";

    TString line;
    while( line.Length() < 40 ) line += " ";
    std::cout << line << std::flush;
    for(size_t s = 0; s < samples.size(); s++) {
      line = samples[s]->name();
      while( line.Length() < 12 ) line += " ";
      std::cout << line << std::flush;
    }
    std::cout << std::endl;

    line = "Total";
    while( line.Length() < 40 ) line += " ";
    std::cout << line << std::flush;
    for(size_t s = 0; s < samples.size(); s++) {
      line = "";
      line += samples[s]->nTotalEvts();
      while( line.Length() < 12 ) line += " ";
      std::cout << line << std::flush;
    }
    std::cout << std::endl;

    for(int c = 0; c < samples[0]->nCuts(); c++) {
      line = samples[0]->printCutLineScreen(c);
      std::cout << line << std::flush;
      for(size_t s = 0; s < samples.size(); s++) {
	line = "";
	line += samples[s]->nPassedEvts(c);
	while( line.Length() < 12 ) line += " ";
	std::cout << line << std::flush;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void writeCutFlowLaTeX(const std::vector<Sample*> samples, const TString &prefix) {
  if( samples.size() > 0 ) {
    TString name = prefix;
    name += "_CutFlowLaTeX";
    std::cout << "Writing cut flow table to LaTeX file '" << name << "'... " << std::flush;
    ofstream file(name);

    file << "\\begin{center}\n";
    file << "\\begin{tabular}[h]{rl";
    for(size_t s = 0; s < samples.size(); s++) {
      file << "r";
    }
    file << "}\n";
    file << " \\\\hline\n\\hline\n";

    file << " & Cut";
    for(size_t s = 0; s < samples.size(); s++) {
      file << " & " << samples[s]->name();
    }
    file << " \\\\\n\\hline\\hline\n";

    file << " & Total";
    for(size_t s = 0; s < samples.size(); s++) {
      file << " & " << samples[s]->nTotalEvts();
    }
    file << " \\\\\n\\hline\n";

    for(int c = 0; c < samples[0]->nCuts(); c++) {
      file << samples[0]->printCutLineLaTeX(c);
      for(size_t s = 0; s < samples.size(); s++) {
	file << " & " << samples[s]->nPassedEvts(c);
      }
      file << " \\\\\n";
    }
    file << "\\hline\\hline\n";
    file << "\\end{tabular}\n";
    file << "\\end{center}\n";
    
    std::cout << "ok" << std::endl;
  }
}

void run(int nMax = -1) {
  int nSamples = 2;
  std::vector<Sample*> samples(nSamples);

  // The data sample
  samples[0] = new Sample("Data","DiJetTree");
  samples[0]->addFile("/afs/naf.desy.de/user/m/mschrode/lustre/data/MinBias-BeamCommissioning09-Dec14thReReco_v1/MinBias-BeamCommissioning09-2360GeV-Dec14thReReco_v1.root");
  samples[0]->setDrawOptions("PE1");

  // The MC sample
  samples[1] = new Sample("MC","DiJetTree");
  samples[1]->addFile("/afs/naf.desy.de/user/m/mschrode/lustre/mc/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1__1.root");
  samples[1]->addFile("/afs/naf.desy.de/user/m/mschrode/lustre/mc/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1__2.root");
  samples[1]->addFile("/afs/naf.desy.de/user/m/mschrode/lustre/mc/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1__3.root");
  samples[1]->addFile("/afs/naf.desy.de/user/m/mschrode/lustre/mc/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1__4.root");
  samples[1]->addFile("/afs/naf.desy.de/user/m/mschrode/lustre/mc/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1/MinBias-Summer09-DESIGN_3X_V8A_2360GeV-v1__5.root");

  for(int s = 0; s < nSamples; s++) {
    samples[s]->readData(nMax);
    samples[s]->fillDistributions();
  }

  draw(samples,"2360 GeV MinBias","plots/2360GeV");
  printCutFlow(samples);
  writeCutFlowLaTeX(samples,"plots/2360GeV");
}

void firstDataAnalysis() {
  run();
}
