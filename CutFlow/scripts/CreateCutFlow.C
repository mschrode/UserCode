#include <iomanip>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"


class CutFlow {
public:
  CutFlow(const TString &fileName);
  ~CutFlow();

  void print() const;

  void addCut(const TString &moduleName, const TString &description);


private:
  class Cut {
  public:
    Cut(const TString &fileName, const TString &moduleName, const TString &description);
    
    double nPassed() const { return nPassed_; }
    TString description() const { return description_; }
    
  private:
    const TString moduleName_;
    const TString description_;
    
    double nPassed_;
    
    TString histName() const { return moduleName_+"/NumEvts"; }
  };

  const TString fileName_;

  std::vector<Cut*> cuts_;

  double efficiencySimple(double nTotal, double nPassed) const { return nTotal > 0. ? nPassed/nTotal : 0.; }
  double efficiencyBayes(double nTotal, double nPassed) const;
};


CutFlow::Cut::Cut(const TString &fileName, const TString &moduleName, const TString &description) 
  : moduleName_(moduleName), description_(description) {
  
  // Get histogram of passed events
  TFile file(fileName,"READ");
  TH1 *h = 0;
  file.GetObject(histName(),h);
  if( h ) {
    nPassed_ = h->GetEntries();
  } else {
    std::cerr << "ERROR in Cut::Cut: Histogram with name '" << histName() << "' does not exist in file '" << fileName << "'\n.";
  }
  file.Close();
}


CutFlow::CutFlow(const TString &fileName)
  : fileName_(fileName) {};


CutFlow::~CutFlow() {
  for(std::vector<Cut*>::iterator it = cuts_.begin(); it != cuts_.end(); ++it) {
    delete *it;
  }
  cuts_.clear();
}


void CutFlow::addCut(const TString &moduleName, const TString &description) {
  cuts_.push_back(new Cut(fileName_,moduleName,description));
}


double CutFlow::efficiencyBayes(double nTotal, double nPassed) const {
  int intTotal = static_cast<int>(nTotal);
  int intPassed = static_cast<int>(nPassed);
  if( intTotal < nTotal || intPassed < nPassed ) {
    std::cerr << "ERROR in CutFlow::efficiencyBayes(): event numbers are non-integers" << std::endl;
    return 0.;
  }

  TH1* hTotal = new TH1D("hTotal","",1,0,1);
  for(int i = 0; i < intTotal; ++i) {
    hTotal->Fill(0);
  }
  TH1* hPassed = new TH1D("hPassed","",1,0,1);
  for(int i = 0; i < intPassed; ++i) {
    hPassed->Fill(0);
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(1);
  g->BayesDivide(hPassed,hTotal,"w");
  
  return g->GetY()[0];
}



void CutFlow::print() const {
  int width = 15;
  for(std::vector<Cut*>::const_iterator it = cuts_.begin(); it != cuts_.end(); ++it) {
    while( (*it)->description().Length() > width ) width += 2;
  }
  std::cout << "\n\n***********************************************\n";
  for(std::vector<Cut*>::const_iterator it = cuts_.begin(); it != cuts_.end(); ++it) {
    std::cout << std::setw(width) << (*it)->description() << " : " <<  (*it)->nPassed() << std::endl;
  }
  std::cout << "***********************************************\n\n\n";


//   double efficiencySimple(double nTotal, double nPassed) const { nTotal > 0. ? nPassed/nTotal : 0. }
//   double efficiencyBayes(double nTotal, double nPassed) const;

}




void CreateCutFlow() {

  CutFlow* cf = new CutFlow("EventCounters.root");
  cf->addCut("EvtCntTotal","Total");
  cf->addCut("EvtCntHLT","HLT_DiJetAve trigger");
  cf->addCut("EvtCntTPBE","TP+BE filter");
  cf->addCut("EvtCntOneGoodVertex","One good vertex filter");
  cf->addCut("EvtCntNoScraping","No-scraping filter");
  cf->addCut("EvtCntHBHENoise","HBHE noise filter");
  cf->addCut("EvtCntEENoise","EE noise filter");
  cf->addCut("EvtCntBeamHalo","Beam halo filter");
  cf->addCut("EvtCntTrackingFailure","Tracking failure filter");
  cf->addCut("EvtCntInconsistentPFMuon","Inconsistent PF muon filter");
  cf->addCut("EvtCntBadPFMuon","Bad PF muon filter");
  cf->addCut("EvtCntGreedyPFMuon","Greedy PF muon filter");
  cf->addCut("EvtCntPFMuon","PF muon filter");
  cf->addCut("EvtCntPFElectron","PF electron filter");

  cf->print();
}



