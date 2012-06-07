// $Id: SystematicUncertainty.h,v 1.2 2011/03/01 16:52:41 mschrode Exp $

#ifndef SYSTEMATIC_UNCERTAINTY_H
#define SYSTEMATIC_UNCERTAINTY_H

#include <vector>

#include "TGraphAsymmErrors.h"
#include "TString.h"

#include "FitResult.h"
#include "Sample.h"

namespace resolutionFit {

  class SystematicUncertainty;
  
  typedef std::vector<SystematicUncertainty*> SystUncerts;
  typedef std::vector<SystematicUncertainty*>::const_iterator SystUncertIt;

  // -------------------------------------------------------------------------------------
  class SystematicUncertainty {
    
  public:
    SystematicUncertainty(const TString &label, int color, const SampleLabel &nominalSample, FitResult::Type type);
    ~SystematicUncertainty();
    
    TString label() const { return label_; }
    FitResult::Type fitResultType() const { return fitResultType_; }
    SampleLabel nominalSampleLabel() const { return sampleLabel_; }
    int color() const { return color_; }

    unsigned int nSteps() const { return ptMean_.size(); }
    double relUncertDown(unsigned int i) const { return relUncertDown_.at(i); }
    double relUncertUp(unsigned int i) const { return relUncertUp_.at(i); }
    TGraphAsymmErrors* relUncertSteps() const;

    bool isCombined() const { return components_.size(); }
    bool hasComponent(const TString &label) const;
    const SystematicUncertainty* component(const TString &label) const;
    unsigned int nComponents() const { return components_.size(); }
    SystUncertIt componentsBegin() const { return components_.begin(); }
    SystUncertIt componentsEnd() const { return components_.end(); }

    void addComponent(const TString &label, int color, const SampleLabel &nominalSample, FitResult::Type type, const std::vector<double> &ptMean, const std::vector<double> &ptMin, const std::vector<double> &ptMax, const std::vector<double> &relUncertsDown, const std::vector<double> &relUncertsUp);

  private:
    const TString label_;
    const int color_;
    const SampleLabel sampleLabel_;
    const FitResult::Type fitResultType_;

    std::vector<double> ptMean_;
    std::vector<double> ptMin_;
    std::vector<double> ptMax_;
    std::vector<double> relUncertDown_;
    std::vector<double> relUncertUp_;

    SystUncerts components_;

    SystematicUncertainty(const TString &label, int color, const SampleLabel &nominalSample, FitResult::Type type, const std::vector<double> &ptMean, const std::vector<double> &ptMin, const std::vector<double> &ptMax, const std::vector<double> &relUncertsDown, const std::vector<double> &relUncertsUp);
  };
}
#endif
