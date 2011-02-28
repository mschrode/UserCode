// $Id: $

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

    TGraphAsymmErrors* relUncertSteps() const;

    bool isCombined() const { return components_.size(); }
    bool hasComponent(const TString &label) const;
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
