// $Id: Parameters.cc,v 1.23 2010/11/11 12:57:04 mschrode Exp $

#include "Parameters.h"

#include <cassert>
#include <iostream>

#include "TRandom3.h"
#include "TStyle.h"

#include "../util/utils.h"


namespace resolutionFit {
  Parameters::Parameters(double etaMin, double etaMax, double deltaPhi12, const TString &fileNameStdSel, const std::vector<double> &ptBinEdges, const TString &outNamePrefix, ResponseFunction::Type type, FitMode fitMode, BinPt binPt, int verbosity)
    : etaMin_(etaMin), etaMax_(etaMax), deltaPhi12_(deltaPhi12),
      styleMode_(gStyle->GetTitle()),
      fitMode_(fitMode),
      binPt_(binPt),
      verbosity_(verbosity),
      nameStdSel_(fileNameStdSel),
      isData_(false) {

    // OutNamePrefix
    TString strEtaMin = util::toTString(10.*etaMin);
    TString strEtaMax = util::toTString(10.*etaMax);
    while( strEtaMin.Length() < 2 ) strEtaMin = strEtaMin+"0";
    while( strEtaMin.Length() < 2 ) strEtaMax = strEtaMax+"0";
    outNamePrefix_ = outNamePrefix+"Eta"+strEtaMin+"-"+strEtaMax+"_";

    init(ptBinEdges, type);
    print();
  }

  void Parameters::init(const std::vector<double> &ptBinEdges, ResponseFunction::Type type) {
    respFunc_ = new ResponseFunction(type);

    hasTruthSpectra_ = false;
    fitPtGenAsym_ = false;
    hasCorrPtGenAsym_ = false;
    ptBinEdges_ = ptBinEdges;

    nameMCStat_ = "";
    nameMCClosure_ = "";

    trueResPar_ = std::vector<double>(3,0.);
    rand_ = new TRandom3(0);
    
    fitExtrapolatedSigma_ = false;
    fitRatio_ = false;
    startResOffset_ = -1.;
    extendedLegend_ = false;
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


  void Parameters::addPt3Threshold(Pt3Var pt3Variable, double pt3RelMax, const TString &fileName, const TString &spectrumName) {
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
      namesCutVars_.push_back(fileName);
      if( fileName == fileNameStdSel() ) stdSelIdx_ = namesCutVars_.size()-1;
      if( hasTruthSpectra() ) {
	if( spectrumName == "" ) {
	  std::cerr << "ERROR: File name for truth spectra must be specified" << std::endl;
	  exit(-2);
	}
	namesTruthSpectra_.push_back(spectrumName);
      }
    }
  }


  void Parameters::addPt3Bin(Pt3Var pt3Variable, double pt3RelMin, double pt3RelMax, double pt3RelMean, const TString &fileName, const TString &spectrumName) {
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
      namesCutVars_.push_back(fileName);
      if( fileName == fileNameStdSel() ) stdSelIdx_ = namesCutVars_.size()-1;
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


  void Parameters::addSystUncert(const TString &label, const TString &fileNameDown, const TString &fileNameUp) {
    labelSyst_.push_back(label);
    namesSystDown_.push_back(fileNameDown);
    namesSystUp_.push_back(fileNameUp);
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
      if( maxLikeFit ) label = "#sigma / <p^{true}_{T}>";
      else label = "#sqrt{2} #sigma_{A} / p^{ave}_{T}";
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
    if( lumi() < 0 ) label = "MC Statistic";
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
    else if( pt3Var() == Pp3Rel ) label = "p_{||,3}";
    label += " < "+util::toTString(pt3Max(pt3Bin));
    if( pt3Var() == Pt3Abs ) label += " GeV";
    else label += (fitMode()==FitModeMaxLikeFull) ? " <p^{true}_{T}>" : " <p^{ave}_{T}>";
    return label;
  }

  
  TString Parameters::labelPtSoftCut() const {
    TString label = (fitMode()==FitModeMaxLikeFull) ? "p_{||,Soft} < 0.015 <p^{true}_{T}>" : "p_{||,Soft} < 0.015 <p^{ave}_{T}>";
    return label;
  }


  TString Parameters::labelMeas() const {
    TString label = "calo_{}";
    return label;
  }


  TString Parameters::labelTruth() const {
    TString label = "true";
    return label;
  }


  TString Parameters::labelPtMeas() const {
    return "p^{"+labelMeas()+"}_{T}";
  }


  TString Parameters::labelPtGen() const {
    return "p^{gen}_{T}";
  }


  TString Parameters::labelPtTrue() const {
    return "p^{"+labelTruth()+"}_{T}";
  }


  TString Parameters::labelPtRef(const TString &plot) const {
    TString label = "p^{ref}_{T}";
    if( binPt() == BinPtGen ) {
      label = labelPtGen();
    } else if(binPt() == BinPtAve ) {
      if( fitMode() == FitModeMaxLikeSimple ) label = "p^{ave}_{T}";
      else if( fitMode() == FitModeMaxLikeFull ) {
	if( plot == "MaxLike" ) label = labelPtGen();
	else if( plot == "PtAsym" ) label = "p^{ave}_{T}";
      }
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
