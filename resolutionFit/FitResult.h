// $Id: $

#ifndef FIT_RESULT_H
#define FIT_RESULT_H

#include <set>
#include <vector>

#include "TF1.h"
#include "TString.h"

#include "Extrapolation.h"
#include "Measurement.h"


namespace resolutionFit {

  // -------------------------------------------------------------------------------------
  class FitResult {
  public:
    enum Type { FullMaxLike, SimpleMaxLike, PtAsym, PtGenAsym };

    static bool validType(Type type);
    static FitResult* createFitResult(Type type, const std::vector<Measurement*> &meas, unsigned int verbosity);
    static TString toString(FitResult::Type type);

    FitResult(const std::vector<Measurement*> &meas, unsigned int verbosity);
    virtual ~FitResult();

    virtual FitResult::Type fitResultType() const = 0;
    
    double meanPt() const { return meanPt_; }
    double meanPtUncert() const { return meanPtUncert_; }
    //    unsigned int nPtSoftBins() const { return meas_.size(); }
    double ptSoft(unsigned int ptSoftBin) const { return ptSoft_.at(ptSoftBin); }
    double value(unsigned int ptSoftBin) const { return values_.at(ptSoftBin); }
    double statUncert(unsigned int ptSoftBin) const { return statUncerts_.at(ptSoftBin); }
    TF1* extrapolationFunction(const TString &name) const {
      return static_cast<TF1*>(extrapolation_->Clone(name)); }
    double extrapolatedValue() const { return extrapolation_->GetParameter(0); }
    double extrapolatedStatUncert() const { return extrapolation_->GetParError(0); }

    void init();

    
  protected:
    const std::vector<Measurement*> meas_;
    const unsigned int verbosity_;

    double meanPt_;
    double meanPtUncert_;
    std::vector<double> values_;
    std::vector<double> statUncerts_;
    std::vector<double> ptSoft_;
    double extrapolatedValue_;
    double extrapolatedStatUncert_;

    virtual bool initResult() = 0;


  private:
    TF1* extrapolation_;

    bool extrapolate();
  };

  typedef std::set<FitResult::Type> FitResultTypes;
  typedef std::set<FitResult::Type>::const_iterator FitResultTypeIt;


  

  // -------------------------------------------------------------------------------------
  class FitResultFullMaxLike : public FitResult {
  public:
    FitResultFullMaxLike(const std::vector<Measurement*> meas, unsigned int verbosity);

    virtual FitResult::Type fitResultType() const { return FitResult::FullMaxLike; }


  protected:
    virtual bool initResult();
  };
}
#endif
