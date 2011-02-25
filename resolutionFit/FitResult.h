// $Id: FitResult.h,v 1.3 2011/02/21 18:25:46 mschrode Exp $

#ifndef FIT_RESULT_H
#define FIT_RESULT_H

#include <map>
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
    enum Type { FullMaxLikeRel, FullMaxLikeAbs, SimpleMaxLike, PtAsym, PtGenAsym };

    static bool validType(Type type);
    static FitResult* createFitResult(Type type, const std::vector<Measurement*> &meas, unsigned int verbosity);
    static TString toString(FitResult::Type type);

    FitResult(const std::vector<Measurement*> &meas, unsigned int verbosity);
    virtual ~FitResult();

    virtual FitResult::Type fitResultType() const = 0;
    
    double meanPt() const { return meanPt_; }
    double meanPtUncert() const { return meanPtUncert_; }
    double ptSoft(unsigned int ptSoftBin) const { return ptSoft_.at(ptSoftBin); }
    double value(unsigned int ptSoftBin) const { return values_.at(ptSoftBin); }
    double statUncert(unsigned int ptSoftBin) const { return statUncerts_.at(ptSoftBin); }
    TF1* extrapolationFunction(const TString &name) const {
      return static_cast<TF1*>(extrapolation_->Clone(name)); }
    double extrapolatedValue() const { return extrapolatedValue_; }
    double extrapolatedStatUncert() const { return extrapolatedStatUncert_; }
    double extrapolatedSystUncert() const { return extrapolatedSystUncert_; }

    virtual bool init() = 0;

    
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
    double extrapolatedSystUncert_;
    TF1* extrapolation_;
  };

  typedef std::set<FitResult::Type> FitResultTypes;
  typedef std::set<FitResult::Type>::const_iterator FitResultTypeIt;


  

  // -------------------------------------------------------------------------------------
  class FitResultFullMaxLikeRel : public FitResult {
  public:
    FitResultFullMaxLikeRel(const std::vector<Measurement*> meas, unsigned int verbosity);

    FitResult::Type fitResultType() const { return FitResult::FullMaxLikeRel; }

    bool init();

  private:
    bool extrapolate();
  };



  // -------------------------------------------------------------------------------------
  class FitResultFullMaxLikeAbs : public FitResult {
  public:
    FitResultFullMaxLikeAbs(const std::vector<Measurement*> meas, unsigned int verbosity);

    FitResult::Type fitResultType() const { return FitResult::FullMaxLikeAbs; }

    bool init();
  };



  // -------------------------------------------------------------------------------------
  class FitResultPtAsym : public FitResult {
  public:
    FitResultPtAsym(const std::vector<Measurement*> meas, unsigned int verbosity);

    FitResult::Type fitResultType() const { return FitResult::PtAsym; }

    bool init();

  private:
    bool extrapolate();
  };
}
#endif
