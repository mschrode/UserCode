// $Id: Parameters.cc,v 1.2 2010/05/18 12:05:43 mschrode Exp $

#include "Parameters.h"

#include <cassert>

#include "TRandom3.h"


namespace resolutionFit {
  Parameters::Parameters(double etaMin, double etaMax, const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, int startIdx, int endIdx, const TString &outNamePrefix, Parameters::ResponseFunction respFunc, int verbosity)
    : etaMin_(etaMin), etaMax_(etaMax),
      startIdx_(startIdx), endIdx_(endIdx),
      outNamePrefix_(outNamePrefix), respFunc_(respFunc), verbosity_(verbosity) {

    ptBinEdges_ = ptBinEdges;
    trueResPar_ = std::vector<double>(3,0.);
    rand_ = new TRandom3(0);

    fitExtrapolatedSigma_ = false;
    fitRatio_ = false;
    startResOffset_ = -1.;

    writeFileNames(namesStdSel_,fileBaseNameStdSel);

    assert( nPtBins() > 0 );
    assert( (endIdx-startIdx) == nPtBins()-1 );
  }

  Parameters::~Parameters() {
    for(std::list<PtBinParameters*>::iterator it = listOfPtBinParameters_.begin();
	it != listOfPtBinParameters_.end(); ++it) {
      delete *it;
    }
    delete rand_;
  }


  const Parameters::PtBinParameters* Parameters::createPtBinParameters(int ptBinIdx) const {
    PtBinParameters *par = new PtBinParameters(ptBinIdx,this);
    listOfPtBinParameters_.push_back(par);
    return par;
  }


  void Parameters::addPt3Cut(double pt3RelCutValue, const TString &fileBaseName) {
    pt3RelCutValues_.push_back(pt3RelCutValue);
    std::vector<TString> names;
    writeFileNames(names,fileBaseName);
    namesCutVars_.push_back(names);
  }


  void Parameters::addSystUncert(const TString &label, const TString &fileBaseNameDown, const TString &fileBaseNameUp) {
    labelSyst_.push_back(label);
    std::vector<TString> names;
    writeFileNames(names,fileBaseNameDown);
    namesSystDown_.push_back(names);
    writeFileNames(names,fileBaseNameUp);
    namesSystUp_.push_back(names);
  }


  void Parameters::addMCTruthBins(int nBins, double min, double max, double relUncert) {
    mcTruthRelUncert_ = relUncert;
    double dPt = (max-min)/nBins;
    for(int i = 0; i <= nBins; ++i) {
      mcTruthPtBinEdges_.push_back(min+i*dPt);
    }
    for(int i = 0; i < nBins; ++i) {
      mcTruthPseudoMeas_.push_back(rand_->Gaus(1.,relUncert));
    }
  }


  void Parameters::writeFileNames(std::vector<TString> &names, const TString &baseName) const {
    names.clear();
    for(int bin = startIdx_; bin <= endIdx_; ++bin) {
      names.push_back(baseName);
      names.back() += bin;
      names.back() += "/jsResponse.root";
    }
  }


  int Parameters::nFittedPars() const {
    int n = 0;
    if( respFunc() == Gauss ) n = 1;
    else if( respFunc() == CrystalBall ) n = 3;

    return n;
  }


  bool Parameters::isRelParValue(int parIdx) const {
    assert( parIdx >=0 && parIdx < nFittedPars() );
    
    bool isRelVal = true;
    if( respFunc() == CrystalBall && parIdx > 0 ) isRelVal = false;

    return isRelVal;
  }

  
  TString Parameters::parLabel(int parIdx) const {
    assert( parIdx >=0 && parIdx < nFittedPars() );

    TString label = "";
    if( parIdx == 0 ) label = "#sigma / p_{T}";
    else if( parIdx == 1 ) label = "#alpha";
    else if( parIdx == 2 ) label = "n";

    return label;
  }


  TString Parameters::parAxisLabel(int parIdx) const {
    assert( parIdx >=0 && parIdx < nFittedPars() );

    TString label = "";
    if( parIdx == 0 ) label = "#sigma / p_{T}";
    else if( parIdx == 1 ) label = "Start  of  tail  #alpha";
    else if( parIdx == 2 ) label = "Slope  of  tail  n";

    return label;
  }


  TString Parameters::labelEtaBin() const {
    char label[50];
    if( etaMin() == 0. ) sprintf(label,"|#eta| < %.1f",etaMax());
    else sprintf(label,"%.1f < |#eta| < %.1f",etaMin(),etaMax());
    
    return label;
  }


  TString Parameters::labelJetAlgo() const {
    return "anti-k_{T} d = 0.5";
  }

  
  TString Parameters::labelPtBin(int ptBin, int type) const {
    char label[50];
    if( type == 0 )      sprintf(label,"%.0f < p^{reco}_{T} < %.0f GeV",ptMin(ptBin),ptMax(ptBin));
    else if( type == 1 ) sprintf(label,"%.0f < p^{particle}_{T} < %.0f GeV",ptMin(ptBin),ptMax(ptBin));
    
    return label;
  }
}
