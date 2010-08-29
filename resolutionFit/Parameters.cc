// $Id: Parameters.cc,v 1.18 2010/08/29 15:54:59 mschrode Exp $

#include "Parameters.h"

#include <cassert>
#include <iostream>

#include "TRandom3.h"
#include "TStyle.h"

#include "../util/utils.h"


namespace resolutionFit {
  Parameters::Parameters(double etaMin, double etaMax, double deltaPhi12, const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, int startIdx, int endIdx, const TString &outNamePrefix, ResponseFunction::Type type, FitMode fitMode, BinPt binPt, int verbosity)
    : etaMin_(etaMin), etaMax_(etaMax), deltaPhi12_(deltaPhi12),
      outNamePrefix_(outNamePrefix),
      styleMode_(gStyle->GetTitle()),
      fitMode_(fitMode),
      binPt_(binPt),
      verbosity_(verbosity),
      isData_(false) {

    for(int i = startIdx; i <= endIdx; ++i) {
      fileNameIdx_.push_back(i);
    }
    init(fileBaseNameStdSel, ptBinEdges, type);

    assert( static_cast<int>(fileNameIdx_.size()) == nPtBins() );
    print();
  }

  Parameters::Parameters(double etaMin, double etaMax, double deltaPhi12, const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, const std::vector<int> fileNameIdx, const TString &outNamePrefix, ResponseFunction::Type type, FitMode fitMode, BinPt binPt, int verbosity)
    : etaMin_(etaMin), etaMax_(etaMax), deltaPhi12_(deltaPhi12),
      outNamePrefix_(outNamePrefix), 
      styleMode_(gStyle->GetTitle()),
      fitMode_(fitMode),
      binPt_(binPt),
      verbosity_(verbosity),
      isData_(false) {
    
    fileNameIdx_ = fileNameIdx;
    init(fileBaseNameStdSel, ptBinEdges, type);

    assert( static_cast<int>(fileNameIdx_.size()) == nPtBins() );
    print();
  }


  void Parameters::init(const TString &fileBaseNameStdSel, const std::vector<double> &ptBinEdges, ResponseFunction::Type type) {
    respFunc_ = new ResponseFunction(type);

    hasTruthSpectra_ = false;
    fitPtGenAsym_ = false;
    hasCorrPtGenAsym_ = false;
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


  void Parameters::addPt3Threshold(Pt3Var pt3Variable, double pt3RelMax, const TString &fileBaseName, const TString &spectrumName) {
    if( nPt3Cuts() == 0 ) {
      pt3Bins_ = false;
      if( spectrumName == "" ) hasTruthSpectra_ = false;
      else hasTruthSpectra_ = true;
      pt3Var_ = pt3Variable;
    }
    if( pt3Bins() ) {
      std::cerr << "ERROR: Already exclusive pt3 binning chosen" << std::endl;
      exit(-2);
    } else if( pt3Variable != pt3Var() ) {
      std::cerr << "ERROR: Already pt3 variable '" << pt3Var() << "' chosen" << std::endl;
      exit(-2);
    } else {
      pt3Max_.push_back(pt3RelMax);
      std::vector<TString> names;
      writeFileNames(names,fileBaseName);
      namesCutVars_.push_back(names);
      if( names.at(0) == fileNameStdSel(0) ) stdSelIdx_ = namesCutVars_.size()-1;
      if( hasTruthSpectra() ) {
	if( spectrumName == "" ) {
	  std::cerr << "ERROR: File name for truth spectra must be specified" << std::endl;
	  exit(-2);
	}
	namesTruthSpectra_.push_back(spectrumName);
      }
    }
  }


  void Parameters::addPt3Bin(Pt3Var pt3Variable, double pt3RelMin, double pt3RelMax, double pt3RelMean, const TString &fileBaseName, const TString &spectrumName) {
    if( nPt3Cuts() == 0 ) {
      pt3Bins_ = true;
      if( spectrumName == "" ) hasTruthSpectra_ = false;
      else hasTruthSpectra_ = true;
      pt3Var_ = pt3Variable;
    }
    if( !pt3Bins() ) {
      std::cerr << "ERROR: Already accumulative pt3 binning chosen" << std::endl;
      exit(-2);
    } else if( pt3Variable != pt3Var() ) {
      std::cerr << "ERROR: Already pt3 variable '" << pt3Var() << "' chosen" << std::endl;
      exit(-2);
    } else {
      pt3Min_.push_back(pt3RelMin);
      pt3Max_.push_back(pt3RelMax);
      pt3Mean_.push_back(pt3RelMean);
      std::vector<TString> names;
      writeFileNames(names,fileBaseName);
      namesCutVars_.push_back(names);
      if( names.at(0) == fileNameStdSel(0) ) stdSelIdx_ = namesCutVars_.size()-1;
      if( nPt3Cuts() == 0 ) {
	if( spectrumName == "" ) hasTruthSpectra_ = false;
	else hasTruthSpectra_ = true;
      }
      if( hasTruthSpectra() ) {
	if( spectrumName == "" ) {
	  std::cerr << "ERROR: File name for truth spectra must be specified" << std::endl;
	  exit(-2);
	}
	namesTruthSpectra_.push_back(spectrumName);
      }
    }
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


  void Parameters::fitPtGenAsym(bool fit) {
    fitPtGenAsym_ = fit;
    hasCorrPtGenAsym_ = fit;
    if( fitPtGenAsym() && ptGenAsymPar_.size() > 0 ) std::cerr << "WARNING: Parameters for ptGen asymmetry already set. Will be overriden by fit!" << std::endl;
  }


  void Parameters::setParPtGenAsym(double a0, double a1, double a2) {
    ptGenAsymPar_.push_back(a0);
    ptGenAsymPar_.push_back(a1);
    ptGenAsymPar_.push_back(a2);
    hasCorrPtGenAsym_ = true;

    if( fitPtGenAsym() ) std::cerr << "WARNING: Parameters for ptGen asymmetry will be overriden by fit!" << std::endl;
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

  
  TString Parameters::parLabel(int parIdx, bool maxLikeFit) const {
    assert( parIdx >=0 && parIdx < nFittedPars() );

    TString label = "";
    if( parIdx == 0 ) {
      if( maxLikeFit ) label = "#sigma / p_{T}";
      else label = "#sqrt{2} #sigma_{A}";
    }
    else if( parIdx == 1 ) label = "#alpha";
    else if( parIdx == 2 ) label = "n";

    return label;
  }


  TString Parameters::parLabelTex(int parIdx, bool maxLikeFit) const {
    assert( parIdx >=0 && parIdx < nFittedPars() );

    TString label = "";
    if( parIdx == 0 ) {
      if( maxLikeFit ) label = "\\sigma / \\pt";
      else label = "\\sqrt{2} \\sigma_{A}";
    }
    else if( parIdx == 1 ) label = "\\alpha";
    else if( parIdx == 2 ) label = "n";

    return label;
  }


  TString Parameters::parAxisLabel(int parIdx) const {
    assert( parIdx >=0 && parIdx < nFittedPars() );

    TString label = "";
    if( parIdx == 1 ) label = "Start  of  tail  ";
    else if( parIdx == 2 ) label = "Slope  of  tail  ";
    label += parLabel(parIdx);

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
    TString label;
    if( lumi() < 0 ) label = "MC Stats";
    else label = "L = "+util::toTString(lumi())+" pb^{-1}";
    return label;
  }

  
  TString Parameters::labelPtBin(int ptBin) const {
    TString label = util::toTString(ptMin(ptBin))+" < p^{";
    if( binPt() == BinPtGen ) label += labelTruth();
    else if( binPt() == BinPtAve ) label += "ave";
    label += "}_{T} < "+util::toTString(ptMax(ptBin))+" GeV";

    return label;
  }


  TString Parameters::labelPt3Cut(int pt3Bin) const {
    TString label;
    if( pt3Bins() ) label += util::toTString(pt3Min(pt3Bin))+" < ";
    if( pt3Var() == Pt3Rel ) label += "p^{rel}_{T,3}";
    else if( pt3Var() == Pt3Abs ) label += "p_{T,3}";
    label += " < "+util::toTString(pt3Max(pt3Bin));
    if( pt3Var() == Pt3Abs ) label += " GeV";
    return label;
  }


  TString Parameters::labelMeas() const {
    TString label = "calo_{}";
    return label;
  }


  TString Parameters::labelTruth() const {
    TString label = "particle";
    return label;
  }


  TString Parameters::labelPtMeas() const {
    return "p^{"+labelMeas()+"}_{T}";
  }


  TString Parameters::labelPtGen() const {
    return "p^{"+labelTruth()+"}_{T}";
  }


  TString Parameters::labelPtRef(const TString &plot) const {
    TString label = "p^{ref}_{T}";
    if( fitMode() == FitModeMaxLikeSimple ) label = "p^{ave}_{T}";
    else if( fitMode() == FitModeMaxLikeFull ) {
      if( plot == "MaxLike" ) label = labelPtGen();
      else if( plot == "PtAsym" ) label = "p^{ave}_{T}";
    }
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


  void Parameters::print() const {
    std::cout << "Setup:\n";
    std::cout << "  FitMode : " << (fitMode()==FitModeMaxLikeSimple ? "FitModeMaxLikeSimple" : "FitModeMaxLikeFull") << std::endl;
    std::cout << "  BinPt   : " << (binPt()==BinPtGen ? "PtGen" : "PtAve") << std::endl;
  }
}
