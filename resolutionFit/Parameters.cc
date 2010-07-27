// $Id: Parameters.cc,v 1.9 2010/07/21 10:55:45 mschrode Exp $

#include "Parameters.h"

#include <cassert>

#include "TRandom3.h"
#include "TStyle.h"

#include "../util/utils.h"


namespace resolutionFit {
  Parameters::Parameters(double etaMin, double etaMax, const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, int startIdx, int endIdx, const TString &outNamePrefix, ResponseFunction::Type type, int verbosity)
    : etaMin_(etaMin), etaMax_(etaMax),
      outNamePrefix_(outNamePrefix),
      styleMode_(gStyle->GetTitle()), verbosity_(verbosity) {

    for(int i = startIdx; i <= endIdx; ++i) {
      fileNameIdx_.push_back(i);
    }
    init(fileBaseNameStdSel, ptBinEdges, type);

    assert( static_cast<int>(fileNameIdx_.size()) == nPtBins() );
  }

  Parameters::Parameters(double etaMin, double etaMax, const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, const std::vector<int> fileNameIdx, const TString &outNamePrefix, ResponseFunction::Type type, int verbosity)
    : etaMin_(etaMin), etaMax_(etaMax),
      outNamePrefix_(outNamePrefix), 
      styleMode_(gStyle->GetTitle()), verbosity_(verbosity) {
    
    fileNameIdx_ = fileNameIdx;
    init(fileBaseNameStdSel, ptBinEdges, type);

    assert( static_cast<int>(fileNameIdx_.size()) == nPtBins() );
  }


  void Parameters::init(const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, ResponseFunction::Type type) {
    respFunc_ = new ResponseFunction(type);
    
    ptBinEdges_ = ptBinEdges;
    ptBinEdges_.erase(ptBinEdges_.begin()+fileNameIdx_.size()+1,ptBinEdges_.end());
    trueResPar_ = std::vector<double>(3,0.);
    rand_ = new TRandom3(0);
    
    fitExtrapolatedSigma_ = false;
    fitRatio_ = false;
    startResOffset_ = -1.;
    extendedLegend_ = false;
    
    writeFileNames(namesStdSel_,fileBaseNameStdSel);
  }


  Parameters::~Parameters() {
    for(std::list<PtBinParameters*>::iterator it = listOfPtBinParameters_.begin();
	it != listOfPtBinParameters_.end(); ++it) {
      delete *it;
    }
    delete respFunc_;
    delete rand_;
  }


  const Parameters::PtBinParameters* Parameters::createPtBinParameters(int ptBinIdx) const {
    PtBinParameters *par = new PtBinParameters(ptBinIdx,this);
    listOfPtBinParameters_.push_back(par);
    return par;
  }


  void Parameters::addPt3Cut(double pt3RelCutValue, const TString &fileBaseName, const TString &spectrumName) {
    pt3RelCutValues_.push_back(pt3RelCutValue);
    std::vector<TString> names;
    writeFileNames(names,fileBaseName);
    namesCutVars_.push_back(names);
    namesTruthSpectra_.push_back(spectrumName);
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
    for(std::vector<int>::const_iterator it = fileNameIdx_.begin();
	it != fileNameIdx_.end(); ++it) {
      names.push_back(baseName);
      names.back() += *it;
      names.back() += "/jsResponse.root";
    }
  }


  bool Parameters::isRelParValue(int parIdx) const {
    assert( parIdx >=0 && parIdx < nFittedPars() );
    
    bool isRelVal = true;
    if( respFuncType() == ResponseFunction::CrystalBall && parIdx > 0 ) isRelVal = false;

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


  TString Parameters::parLabelTex(int parIdx) const {
    assert( parIdx >=0 && parIdx < nFittedPars() );

    TString label = "";
    if( parIdx == 0 ) label = "\\sigma / \\pt";
    else if( parIdx == 1 ) label = "\\alpha";
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


  TString Parameters::labelLumi() const {
    return "L = "+util::toTString(lumi())+" pb^{-1}";
  }

  
  TString Parameters::labelPtBin(int ptBin, int type) const {
    char label[50];
    if( type == 0 )      sprintf(label,"%.0f < p^{%s}_{T} < %.0f GeV",
				 ptMin(ptBin),labelMeas().Data(),ptMax(ptBin));
    else if( type == 1 ) sprintf(label,"%.0f < p^{%s}_{T} < %.0f GeV",
				 ptMin(ptBin),labelTruth().Data(),ptMax(ptBin));    
    return label;
  }


  TString Parameters::labelPt3Cut(int ptBin) const {
    char label[50];
    sprintf(label,"p^{rel}_{T,3} < %.2f",pt3CutValue(ptBin));
    return label;
  }


  TString Parameters::labelMeas() const {
    TString label = "calo_{}";
    if( styleMode() == "CMS" ) label = "Jet";
    return label;
  }


  TString Parameters::labelTruth() const {
    TString label = "particle";
    if( styleMode() == "CMS" ) label = "GenJet";
    return label;
  }


  TString Parameters::xAxisTitleResponse() const {
    TString title = "";
    if( extendedLegend() ) title += "Response ";
    title += "R = p^{" + labelMeas() + "}_{T} / p^{" + labelTruth() + "}_{T}";
    if( styleMode() == "CMS" ) title = "p^{" + labelMeas() + "}_{T} / p^{" + labelTruth() + "}_{T}";
    return title;
  }


  TString Parameters::yAxisTitleResponse() const {
    return "1 / N  dN / dR";
  }

}
