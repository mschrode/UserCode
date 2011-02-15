// $Id: $

#include <iostream>

#include "EtaBin.h"


namespace resolutionFit{

  // -------------------------------------------------------------------------------------
  EtaBin::EtaBin(unsigned int etaBin, unsigned int nPtBins, const Parameters* par)
    : par_(par), etaBin_(etaBin) {
    for(unsigned int i = 0; i < nPtBins; ++i) {
      ptBins_.push_back(new PtBin(etaBin_,i,par));
    }
    mcTruthReso_ = new ResolutionFunctionModifiedNSC(0.,1.,0.,1.,0.,0.);
    pli_ = new ResolutionFunctionModifiedNSC(0.,1.,0.,1.,0.,0.);
  }


  // -------------------------------------------------------------------------------------
  EtaBin::~EtaBin() {
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      delete *it;
    }
    delete pli_;
    delete mcTruthReso_;
  }


  // -------------------------------------------------------------------------------------
  TString EtaBin::toString() const {
    TString str = "Eta";
    str += etaBin(); 

    return str;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addDataSample(const TString &label, const TString &baseFileName) {
    bool result = true;
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      bool status = (*it)->addDataSample(label,baseFileName);
      result = result && status;
    }
    sampleTypes_[label] = Sample::Data;

    return result;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addMCSample(const TString &label, const TString &baseFileName) {
    bool result = true;
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      bool status = (*it)->addMCSample(label,baseFileName);
      result = result && status;
    }
    sampleTypes_[label] = Sample::MC;

    return result;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addFitResult(FitResult::Type type) {
    bool result = true;
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      bool status = (*it)->addFitResult(type);
      result = result && status;
    }
    fitResultTypes_.insert(type);

    return result;
  }


  // -------------------------------------------------------------------------------------
  void EtaBin::setPLI(ResolutionFunction* pli) { 
    delete pli_;
    pli_ = pli;
  }


  // -------------------------------------------------------------------------------------
  void EtaBin::setMCTruthResolution(ResolutionFunction* mcTruthReso) { 
    delete mcTruthReso_;
    mcTruthReso_ = mcTruthReso;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::meanPt(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->meanPt(type);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::meanPtStatUncert(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->meanPtStatUncert(type);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::extrapolatedValue(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->extrapolatedValue(type);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::extrapolatedStatUncert(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->extrapolatedStatUncert(type);
	break;
      }
    }

    return result;
  }


  // Should be integral over bin!
  // -------------------------------------------------------------------------------------
  double EtaBin::mcTruthResolution(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const {
    return mcTruthReso_->val(meanPt(label,type,ptBinIdx));
  }


  // Should be integral over bin!
  // -------------------------------------------------------------------------------------
  double EtaBin::pli(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const {
    return pli_->val(meanPt(label,type,ptBinIdx));
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::correctedResolution(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = extrapolatedValue(label,type,ptBinIdx);
    double plim = pli(label,type,ptBinIdx);

    return sqrt( result*result - plim*plim );
  }

  
  // -------------------------------------------------------------------------------------
  double EtaBin::correctedResolutionStatUncert(SampleLabel label, FitResult::Type type, unsigned int ptBinIdx) const {
    return extrapolatedStatUncert(label,type,ptBinIdx);
  }
}
