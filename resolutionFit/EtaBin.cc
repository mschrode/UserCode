// $Id: EtaBin.cc,v 1.4 2011/02/26 17:55:50 mschrode Exp $

#include <algorithm>
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
  bool EtaBin::addMCClosureUncertainty(const SampleLabel &nominalSample, FitResult::Type type, int color) {
    SystematicUncertainty* uncert = findSystematicUncertainty(nominalSample,type);
    
    TString label = "MC Closure";
    if( !(uncert->hasComponent(label)) ) {
      std::vector<double> ptMean;
      std::vector<double> ptMin;
      std::vector<double> ptMax;
      std::vector<double> relUncert;
      for(unsigned int i = 0; i < par_->nPtBins(etaBin()); ++i) {
 	ptMean.push_back(meanPt(nominalSample,type,i));
 	ptMin.push_back(par_->ptMin(etaBin(),i));
 	ptMax.push_back(par_->ptMax(etaBin(),i));
	relUncert.push_back(0.5*std::abs(1.-relativeMCClosure(nominalSample,type,i)));
      }
      uncert->addComponent(label,color,nominalSample,type,ptMean,ptMin,ptMax,relUncert,relUncert);
    }
    
    return true;
  }


  // -------------------------------------------------------------------------------------
  bool EtaBin::addUncertaintyFromVariedSample(const TString &uncertaintyLabel, double fraction, const SampleLabel &nominalSample, FitResult::Type type, const TString &variedSampleDown, const TString &variedSampleUp, int color) {
    SystematicUncertainty* uncert = findSystematicUncertainty(nominalSample,type);
    
    if( !(uncert->hasComponent(uncertaintyLabel)) ) {
      std::vector<double> ptMean;
      std::vector<double> ptMin;
      std::vector<double> ptMax;
      std::vector<double> relUncert;
      for(unsigned int i = 0; i < par_->nPtBins(etaBin()); ++i) {
 	ptMean.push_back(meanPt(nominalSample,type,i));
 	ptMin.push_back(par_->ptMin(etaBin(),i));
 	ptMax.push_back(par_->ptMax(etaBin(),i));
	double res = correctedResolution(nominalSample,type,i);
	double resDown = correctedResolution(variedSampleDown,type,i);
	double resUp = correctedResolution(variedSampleUp,type,i);
	relUncert.push_back(fraction*0.5*(std::abs((resDown-res)/res)+std::abs((resUp-res)/res)));
      }
      uncert->addComponent(uncertaintyLabel,color,nominalSample,type,ptMean,ptMin,ptMax,relUncert,relUncert);
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
  void EtaBin::fitPLI(const TString &label, const TString &baseFileName, ResolutionFunction::Type type) {

    // Add Samples containing ptGen asymmetry histograms
    // They will per default have the correct FitResult::Type PtGenAsym
    for(std::vector<PtBin*>::iterator it = ptBins_.begin(); it != ptBins_.end(); ++it) {
      if( !((*it)->addMCTruthSample(label,baseFileName)) ) {
	std::cerr << "ERROR in EtaBin::fitPLI(): Error adding MCTruth Sample '" << label << "'" << std::endl;
	exit(1);
      }
    }

    // Create graph of extrapolated ptGen asymmetry width
    std::vector<double> pt;
    std::vector<double> ptErr;
    std::vector<double> res;
    std::vector<double> resStatErr;
    for(PtBinIt ptBinIt = ptBinsBegin(); ptBinIt != ptBinsEnd(); ++ptBinIt) {
      const Sample* sample = (*ptBinIt)->mcTruthSample();
      pt.push_back(sample->meanPt(FitResult::PtGenAsym));
      ptErr.push_back(0.);
      res.push_back(sample->extrapolatedValue(FitResult::PtGenAsym));
      resStatErr.push_back(sample->extrapolatedStatUncert(FitResult::PtGenAsym));
    }
    TGraphAsymmErrors* g = new TGraphAsymmErrors(pt.size(),&(pt.front()),&(res.front()),
						 &(ptErr.front()),&(ptErr.front()),
						 &(resStatErr.front()),&(resStatErr.front()));

    // Fit new PLI from extrapolated ptGen asymmetry width
    delete pli_;
    pli_ = fitResolution(g,type);
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


  //! Returns MCTruth / Measurement
  // -------------------------------------------------------------------------------------
  double EtaBin::relativeMCClosure(const SampleLabel &label, FitResult::Type type, unsigned int ptBinIdx) const {
    double meas = correctedResolution(label,type,ptBinIdx);
    double truth = mcTruthResolution(label,type,ptBinIdx);
    double result = 0.;
    if( meas > 0. ) result = truth/meas;

    return result;
  }


  //! Create graph of extrapolated resolution
  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::extrapolatedResolution(const SampleLabel &label, FitResult::Type type) const {
    std::vector<double> pt;
    std::vector<double> ptErr;
    std::vector<double> res;
    std::vector<double> resStatErr;
    for(PtBinIt ptBinIt = ptBinsBegin(); ptBinIt != ptBinsEnd(); ++ptBinIt) {
      unsigned int ptBinIdx = (*ptBinIt)->ptBin();
      pt.push_back(meanPt(label,type,ptBinIdx));
      ptErr.push_back(0.);
      res.push_back(extrapolatedValue(label,type,ptBinIdx));
      resStatErr.push_back(extrapolatedStatUncert(label,type,ptBinIdx));
    }

    return new TGraphAsymmErrors(pt.size(),&(pt.front()),&(res.front()),&(ptErr.front()),&(ptErr.front()),
				 &(resStatErr.front()),&(resStatErr.front()));
  }



  //! Create graph of extrapolated and PLI-corrected resolution
  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* EtaBin::correctedResolution(const SampleLabel &label, FitResult::Type type) const {
    std::vector<double> pt;
    std::vector<double> ptErr;
    std::vector<double> res;
    std::vector<double> resStatErr;
    for(PtBinIt ptBinIt = ptBinsBegin(); ptBinIt != ptBinsEnd(); ++ptBinIt) {
      unsigned int ptBinIdx = (*ptBinIt)->ptBin();
      pt.push_back(meanPt(label,type,ptBinIdx));
      ptErr.push_back(0.);
      res.push_back(correctedResolution(label,type,ptBinIdx));
      resStatErr.push_back(correctedResolutionStatUncert(label,type,ptBinIdx));
    }

    return new TGraphAsymmErrors(pt.size(),&(pt.front()),&(res.front()),&(ptErr.front()),&(ptErr.front()),
				 &(resStatErr.front()),&(resStatErr.front()));
  }


  // -------------------------------------------------------------------------------------
  ResolutionFunction* EtaBin::fitResolution(const TGraphAsymmErrors* g, ResolutionFunction::Type type) const {

    TGraphAsymmErrors* gc = static_cast<TGraphAsymmErrors*>(g->Clone());

    std::vector<double> param;
    // Min and max pt
    param.push_back(*std::min_element(gc->GetX(),gc->GetX()+gc->GetN()));
    param.push_back(*std::max_element(gc->GetX(),gc->GetX()+gc->GetN()));
    // Choose typical start values
    if( type == ResolutionFunction::NSC ) {
      param.push_back(4.);
      param.push_back(1.);
      param.push_back(0.03);
    } else if( type == ResolutionFunction::ModifiedNSC ) {
      param.push_back(4.);
      param.push_back(0.7);
      param.push_back(0.);
      param.push_back(0.2);
    }
    ResolutionFunction* f = ResolutionFunction::createResolutionFunction(type,param);

    // Get TF1 to fit
    TF1* fit = f->func("tmp");
    if( type == ResolutionFunction::ModifiedNSC ) fit->FixParameter(2,0.);
    gc->Fit(fit,"0BR");
    ResolutionFunction* result = ResolutionFunction::createResolutionFunction(type,fit);
    delete f;
    delete fit;
    delete gc;
    
    return result;
  }
}
