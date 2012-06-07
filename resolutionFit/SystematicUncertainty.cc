// $Id: SystematicUncertainty.cc,v 1.2 2012/05/31 20:17:44 mschrode Exp $

#include <cassert>
#include <cmath>

#include "SystematicUncertainty.h"


namespace resolutionFit {

  //! Constructor for combined systematic uncertainty
  //! Use addComponent to add individual components
  // -------------------------------------------------------------------------------------
  SystematicUncertainty::SystematicUncertainty(const TString &label, int color, const SampleLabel &nominalSample, FitResult::Type type) 
    : label_(label), color_(color), sampleLabel_(nominalSample), fitResultType_(type) {};


  //! Constructor for individual component
  // -------------------------------------------------------------------------------------
  SystematicUncertainty::SystematicUncertainty(const TString &label, int color, const SampleLabel &nominalSample, FitResult::Type type, const std::vector<double> &ptMean, const std::vector<double> &ptMin, const std::vector<double> &ptMax, const std::vector<double> &relUncertsDown, const std::vector<double> &relUncertsUp) 
    : label_(label), color_(color), sampleLabel_(nominalSample), fitResultType_(type) {

    assert( ptMean.size() == ptMin.size() );
    assert( ptMean.size() == ptMax.size() );
    assert( ptMean.size() == relUncertsDown.size() );
    assert( ptMean.size() == relUncertsUp.size() );
    
    for(size_t i = 0; i < ptMean.size(); ++i) {
      ptMean_.push_back(ptMean.at(i));
      ptMin_.push_back(ptMin.at(i));
      ptMax_.push_back(ptMax.at(i));
      relUncertDown_.push_back(relUncertsDown.at(i));
      relUncertUp_.push_back(relUncertsUp.at(i));
    }
  }


  // -------------------------------------------------------------------------------------
  SystematicUncertainty::~SystematicUncertainty() {
    for(SystUncerts::iterator it = components_.begin(); it != components_.end(); ++it) {
      delete *it;
    }
  }


  // -------------------------------------------------------------------------------------
  void SystematicUncertainty::addComponent(const TString &label, int color, const SampleLabel &nominalSample, FitResult::Type type, const std::vector<double> &ptMean, const std::vector<double> &ptMin, const std::vector<double> &ptMax, const std::vector<double> &relUncertsDown, const std::vector<double> &relUncertsUp) {

    if( isCombined() ) {
      // Add uncertainty in quadrature to total uncertainty
      for(size_t i = 0; i < relUncertDown_.size(); ++i) {
	relUncertDown_.at(i) = sqrt( relUncertDown_.at(i)*relUncertDown_.at(i) + relUncertsDown.at(i)*relUncertsDown.at(i) );
	relUncertUp_.at(i) = sqrt( relUncertUp_.at(i)*relUncertUp_.at(i) + relUncertsUp.at(i)*relUncertsUp.at(i) );
      }
    } else {
      ptMean_ = ptMean;
      ptMin_ = ptMin;
      ptMax_ = ptMax;
      relUncertDown_ = relUncertsDown;
      relUncertUp_ = relUncertsUp;
    }

    // Add component to list of components
    components_.push_back(new SystematicUncertainty(label,color,nominalSample,type,ptMean,ptMin,ptMax,relUncertsDown,relUncertsUp));
  }


  // -------------------------------------------------------------------------------------
  TGraphAsymmErrors* SystematicUncertainty::relUncertSteps() const {
    std::vector<double> x;
    std::vector<double> xed;
    std::vector<double> xeu;
    std::vector<double> y;
    std::vector<double> yed;
    std::vector<double> yeu;
    for(size_t i = 0; i < ptMean_.size(); ++i) {
      x.push_back(ptMean_.at(i));
      xed.push_back(x.at(i)-ptMin_.at(i));
      xed.push_back(ptMax_.at(i)-x.at(i));
      y.push_back(0.);
      yed.push_back(relUncertDown_.at(i));
      yeu.push_back(relUncertUp_.at(i));
    }
    TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(),&(x.front()),&(y.front()),
						 &(xed.front()),&(xeu.front()),&(yed.front()),&(yeu.front()));
    g->SetFillStyle(1001);
    g->SetFillColor(color_);
    g->SetLineColor(color_);

    return g;
  }


  // -------------------------------------------------------------------------------------
  bool SystematicUncertainty::hasComponent(const TString &label) const {
    bool result = false;
    for(SystUncertIt it = components_.begin(); it != components_.end(); ++it) {
      if( (*it)->label() == label ) {
	result = true;
	break;
      }
    }
    
    return result;
  }


  // -------------------------------------------------------------------------------------
  const SystematicUncertainty* SystematicUncertainty::component(const TString &label) const {
    const SystematicUncertainty* s = 0;
    for(SystUncertIt it = components_.begin(); it != components_.end(); ++it) {
      if( (*it)->label() == label ) {
	s = *it;
	break;
      }
    }
    
    return s;
  }
}
