// $Id: EtaBin.cc,v 1.2 2011/02/18 18:42:22 mschrode Exp $

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
    for(ComparedSamples::iterator it = compSamples_.begin(); it != compSamples_.end(); ++it) {
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
  bool EtaBin::hasSystematicUncertainty(const SampleLabel &label, FitResult::Type type) const {
    bool result = false;
    for(std::set<SystematicUncertainty*>::const_iterator it = systUncerts_.begin();
	it != systUncerts_.end(); ++it) {
      if( (*it)->nominalSampleLabel() == label && (*it)->fitResultType() == type ) {
	result = true;
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::findSystematicUncertainty(const SampleLabel &label, FitResult::Type type, const SystematicUncertainty* &uncert) const {
    uncert = 0;
    bool result = false;
    for(std::set<SystematicUncertainty*>::const_iterator it = systUncerts_.begin();
	it != systUncerts_.end(); ++it) {
      if( (*it)->nominalSampleLabel() == label && (*it)->fitResultType() == type ) {
	result = true;
	uncert = *it;
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  SystematicUncertainty* EtaBin::findSystematicUncertainty(const SampleLabel &label, FitResult::Type type) {
    SystematicUncertainty* uncert = 0;

    // Check whether there are uncertainties already associated to this Sample and FitResult
    for(std::set<SystematicUncertainty*>::const_iterator it = systUncerts_.begin();
	it != systUncerts_.end(); ++it) {
      if( (*it)->nominalSampleLabel() == label && (*it)->fitResultType() == type ) {
	uncert = *it;
	break;
      }
    }

    // If not, create new combined uncertainty and add to others
    if( uncert == 0 ) {
      uncert = new SystematicUncertainty("Total",5,label,type);
      systUncerts_.insert(uncert);
    }

    return uncert;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addExtrapolationUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color) {
    SystematicUncertainty* uncert = findSystematicUncertainty(nominalSample,type);

    TString label = "Extrapolation";
    if( !(uncert->hasComponent(label)) ) {
      std::vector<double> ptMean;
      std::vector<double> ptMin;
      std::vector<double> ptMax;
      std::vector<double> relUncerts;
      for(unsigned int i = 0; i < par_->nPtBins(etaBin()); ++i) {
	ptMean.push_back(meanPt(nominalSample,type,i));
	ptMin.push_back(par_->ptMin(etaBin(),i));
	ptMax.push_back(par_->ptMax(etaBin(),i));
	relUncerts.push_back(std::abs(extrapolatedSystUncert(nominalSample,type,i)/extrapolatedValue(nominalSample,type,i)));
      }
      uncert->addComponent(label,color,nominalSample,type,ptMean,ptMin,ptMax,relUncerts,relUncerts);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addPLIUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color) {
    SystematicUncertainty* uncert = findSystematicUncertainty(nominalSample,type);
    
    TString label = "Particle Level Imbalance";
    if( !(uncert->hasComponent(label)) ) {
      std::vector<double> ptMean;
      std::vector<double> ptMin;
      std::vector<double> ptMax;
      std::vector<double> relUncertDown;
      std::vector<double> relUncertUp;
      for(unsigned int i = 0; i < par_->nPtBins(etaBin()); ++i) {
 	ptMean.push_back(meanPt(nominalSample,type,i));
 	ptMin.push_back(par_->ptMin(etaBin(),i));
 	ptMax.push_back(par_->ptMax(etaBin(),i));
	double res = correctedResolution(nominalSample,type,i);
	relUncertDown.push_back(std::abs((correctedResolution(nominalSample,type,i,0.75)-res)/res));
	relUncertUp.push_back(std::abs((correctedResolution(nominalSample,type,i,1.25)-res)/res));
      }
      uncert->addComponent(label,color,nominalSample,type,ptMean,ptMin,ptMax,relUncertDown,relUncertUp);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::compareSamples(const SampleLabel &label1, const SampleLabel &label2) {
    if( sampleTypes_.find(label1) != sampleTypes_.end() &&
	sampleTypes_.find(label2) != sampleTypes_.end() ) {
      compSamples_.insert(new SampleLabelPair(label1,label2));
    } else {
      std::cerr << "ERROR in EtaBin::compareSamples(): no samples with label '" << label1 << "' and '" << label2 << "'" << std::endl;
      exit(1);
    }

    return true;
  }


  // -------------------------------------------------------------------------------------
  Sample::Type EtaBin::sampleType(const SampleLabel &label) const {
    SampleTypeIt it = sampleTypes_.find(label);
    if( it == sampleTypes_.end() ) {
      std::cerr << "ERROR in EtaBin::sampleType(): No sample with label '" << label << "'" << std::endl;
      exit(1);
    }

    return it->second;
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
  double EtaBin::meanPt(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
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
  double EtaBin::meanPtStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
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
  double EtaBin::extrapolatedValue(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
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
  double EtaBin::extrapolatedStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->extrapolatedStatUncert(type);
	break;
      }
    }

    return result;
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::extrapolatedSystUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    double result = 0.;
    for(PtBinIt ptBinIt = ptBins_.begin(); ptBinIt != ptBins_.end(); ++ptBinIt) {
      if( (*ptBinIt)->ptBin() == ptBinIdx ) {
	result = (*ptBinIt)->findSample(label)->extrapolatedSystUncert(type);
	break;
      }
    }

    return result;
  }


  // Should be integral over bin!
  // -------------------------------------------------------------------------------------
  double EtaBin::mcTruthResolution(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    return mcTruthReso_->val(meanPt(label,type,ptBinIdx));
  }


  // Should be integral over bin!
  // -------------------------------------------------------------------------------------
  double EtaBin::pli(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    return pli_->val(meanPt(label,type,ptBinIdx));
  }


  // -------------------------------------------------------------------------------------
  double EtaBin::correctedResolution(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx, double scalePLI) const {
    double result = extrapolatedValue(label,type,ptBinIdx);
    double plim = scalePLI*pli(label,type,ptBinIdx);

    return sqrt( result*result - plim*plim );
  }

  
  // -------------------------------------------------------------------------------------
  double EtaBin::correctedResolutionStatUncert(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    return extrapolatedStatUncert(label,type,ptBinIdx);
  }
}
