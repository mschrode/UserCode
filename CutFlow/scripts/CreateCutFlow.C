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

  void print(bool latex = false, int prescision = 5) const;

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
  void efficiencyBayes(double nTotal, double nPassed, double &eff, double &effErrDown, double &effErrUp) const;
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


void CutFlow::efficiencyBayes(double nTotal, double nPassed, double &eff, double &effErrDown, double &effErrUp) const {
  int intTotal = static_cast<int>(nTotal);
  int intPassed = static_cast<int>(nPassed);
  if( intTotal < nTotal || intPassed < nPassed ) {

    std::cerr << "ERROR in CutFlow::efficiencyBayes(): event numbers are non-integers" << std::endl;
    eff = 0.;
    effErrDown = 0.;
    effErrUp = 0.;

  } else {

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
    
    eff = g->GetY()[0]; 
    effErrDown = g->GetEYlow()[0];
    effErrUp = g->GetEYhigh()[0];

    delete hTotal;
    delete hPassed;
    delete g;
  }
}



void CutFlow::print(bool latex, int prescision) const {
  int width = 15;
  for(std::vector<Cut*>::const_iterator it = cuts_.begin(); it != cuts_.end(); ++it) {
    while( (*it)->description().Length() > width ) width += 2;
  }
  if( latex ) {
    std::cout << "\n\n\\begin{tabular}{rrrr}\n\\hline\n";
    std::cout << std::setw(width) << "Cut     ";
    std::cout << " & ";
    std::cout << std::setw(10) << "$N(\\text{passed})$";
    std::cout << " & ";
    std::cout << std::setw(10) << "$\epsilon(\\text{binomial})$";
    std::cout << " & ";
    std::cout << std::setw(34) << "$\epsilon(\\text{bayesian})$ \\\\ \n";
  } else {
    std::cout << "\n\n**********************************************************************************************\n";
    std::cout << std::setw(width) << "Cut     ";
    std::cout << " : ";
    std::cout << std::setw(10) << "# passed ";
    std::cout << " : ";
    std::cout << std::setw(10) << "Eff   ";
    std::cout << " : ";
    std::cout << std::setw(34) << "Eff Bayes                \n";
    std::cout << "----------------------------------------------------------------------------------------------\n";
  }
  for(size_t i = 0; i < cuts_.size(); ++i) {
    if( latex ) {
      std::cout << std::setw(width) << "\\texttt{" << cuts_.at(i)->description() << "} & $";
      std::cout << std::setw(10) << cuts_.at(i)->nPassed() << "$ & $";
    } else {
      std::cout << std::setw(width) << cuts_.at(i)->description() << " : ";
      std::cout << std::setw(10) << cuts_.at(i)->nPassed() << " : ";
    }
    if( i > 0 ) {
      double nTotal = cuts_.at(i-1)->nPassed();
      double nPassed = cuts_.at(i)->nPassed();
      
      double eff = 0.;
      double effErrDown = 0.;
      double effErrUp = 0.;
      //efficiencyBayes(nTotal,nPassed,eff,effErrDown,effErrUp);

      if( latex ) {
	std::cout.precision(prescision);
	std::cout << std::setw(10) << std::fixed <<  efficiencySimple(nTotal,nPassed) << "$ & $";
	std::cout << std::setw(10) << eff << "^{+" << effErrUp << "_{-" << effErrDown << "}";
      } else {
	std::cout << std::setw(10) << efficiencySimple(nTotal,nPassed) << " : ";
	std::cout << std::setw(10) << eff << " -" << effErrDown << " +" << effErrUp;
      }
    }
    if( latex ) std::cout << " \\\\";
    std::cout << std::endl;
  }
  if( !latex) std::cout << "**********************************************************************************************";
  std::cout << "\n\n\n";
}




void CreateCutFlow(const TString &fileName, bool latex = false, int precision = 5) {

  CutFlow* cf = new CutFlow(fileName);
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

  cf->print(latex,precision);
}



